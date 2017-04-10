/*
*  Author : Jeffrey A. Robinson
*  Date   : August 20, 2013
*
*  Implementation of Systolic Montgomery Multiplication, Addition, and Subtraction.
*  Systolic Montgomery Multiplication : B. Dixon, A.K. Lenstra 
*
*/

#define DEFINES 1
#define GPU_HELPER_FUNCTIONS 1
#define CPU_HELPER_FUNCTIONS 1
#define EXTERNS 1
#define PROTOTYPES 1

// Stdlib
#include <string>
#include <iostream>
#include <random>

// My header files and stuff
#include "../Integer/Integer.h"

// Cuda
#include <cuda_runtime.h>
#include <cuda_texture_types.h>

// Use these namespaces.
using namespace std;

extern const uint3 threadIdx;
extern const uint3 blockIdx;
extern const dim3 gridDim;
extern const dim3 blockDim;

#if DEFINES
#define _VOL_
#define IDX threadIdx.x
#define BDX blockIdx.x
#endif

#if EXTERNS
// Make intellisense leave me alone

// atomics and syncing
extern __device__ void __syncthreads();
extern __device__ int __syncthreads_count(int predicate);
extern __device__ int __syncthreads_or(int predicate);
extern __device__ int __syncthreads_and(int predicate);
extern __device__ void __threadfence();
extern __device__  int atomicMax(  int* address ,  int val );
extern __device__  int atomicMin(  int* address ,  int val );
extern __device__  unsigned int atomicOr(  unsigned int* address ,  unsigned int val );

// Bit Conversions
extern __device__ double __longlong_as_double( long long x );
extern __device__ int __float_as_int( int x );

extern __device__ long long __double_as_longlong( double x );
extern __device__ float __int_as_float( float x );

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		Do not use this, it is possible where your neighbors will not get the shifts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  */
// New operation in CUDA 3.0+.  Just like MASPAR operations where you can shift 
// to processors to the "east" and "west"
extern __device__ int __shfl_up( int , unsigned int , unsigned int );
extern __device__ int __shfl_down( int , unsigned int , unsigned int );
#endif


#if PROTOTYPES
__device__ void d_mul( _VOL_ C_P_BASE , _VOL_ C_P_BASE , C_P_BASE , BASE , BASE  );
void validate( Integer*A , Integer*B, Integer *M , Integer*R , int r , int d , Integer* Results , std::ostream );
#endif

#if CPU_HELPER_FUNCTIONS
// print out the error.
#define errorTest(  err , out ) { \
	if( err != cudaSuccess ) { \
		out << cudaGetErrorString( err ) << endl; \
	} \
}

// make sure to only print out errors and not the successes
#define errorTestMsg(  msg , err , out ) { \
	if( err != cudaSuccess ) { \
		out << msg << " : " << cudaGetErrorString( err ) << endl; \
	} \
}
#endif

#if 0
// We only need the bottom 23 bits
#define lowF( f ) __longlong_as_double( __double_as_longlong( f ) & (0x7FFFFF) )

// We only need the bottom 23 bits
#define hiF( f ) __longlong_as_double( (__double_as_longlong( f ) >> 23) & (0x7FFFFF) )
#endif


// Shared Conditional Variable
extern __shared__ BASE TMP[];


#if __CUDA_ARCH__ >= 500
#define shift_up( in ) in = (__shfl_up( int(in & MASK) , (unsigned int)1 , SIZE ) * (IDX != 0))
#define shift_down( in ) in = (__shfl_down( int(in & MASK) , (unsigned int)1 , SIZE ) * (IDX != SIZE))
#elif 0
#define shift_up(su_move,su_overflow) \
{ \
	TMP[(IDX + 1) % r] = low(su_move); \
	__syncthreads(); \
	su_overflow = TMP[0]; \
	su_move = TMP[IDX] * (IDX>0); \
}

#define shift_down( in , shift_in ) { \
	TMP[IDX] = in; \
	TMP[0] = shift_in; \
	__syncthreads(); \
	in = TMP[(IDX + 1)%r]; \
}

#endif

#if 0

/* Helper functions */

#define single_prop(sp_v,sp_o) \
{ \
	TMP[(IDX + 1)%blockDim.x] = hi(sp_v); \
	__syncthreads(); \
	sp_v = low(sp_v); \
	sp_v += (TMP[IDX]) * (IDX != 0); \
	sp_o = TMP[0]; \
	__syncthreads(); \
}

#define broadCast(_v_,id,_r_) { \
	__shared__ unsigned int t[2]; \
	t[ (IDX != id)] = _v_; \
	__syncthreads(); \
	_r_ = t[0]; \
}

// Propogate the carries
#define prop(pr_v,pr_o) \
{ \
	pr_o = 0;\
	do { \
		cond = 0; \
		BASE pr_tmp = 0;\
		single_prop( pr_v , pr_tmp ); \
		pr_o += pr_tmp; \
		if( hi(pr_v) ) cond = 1;\
		__syncthreads(); \
	} while( cond ); \
}

#endif


extern __device__ float	__fdiv_rz( float , float );
extern __device__ int	__float2int_rz( float );
extern __device__ float	__int2float_rz( int );
#define mod(x,y) (x)%(y)
//(x) - (y) * int(float(x)/float(y))


// Works!
__device__ void d_mul_concurrent_sync( P_BASE __restrict A , P_BASE __restrict B , P_BASE __restrict R , BASE n, unsigned int r , unsigned int d ) {
	
	// Used to test
	__shared__ int cond[33];

	// Result will be stored here
	unsigned long long u = 0;

	// Cache locally
	unsigned long long a = A[IDX];
	BASE b = B[IDX];

	int index = IDX / r;
	int base = index * r;
	int up = mod(IDX + 1,r) + base;
	__syncthreads();

	// Perform operations
	for( unsigned int i = 0 ; i < r ; i++ ) {		
		TMP[IDX] = b;
		__syncthreads();
		// Stall here
		u = TMP[base + i] * a + u;	// Do the computation
		TMP[IDX] = u; // save the carry
		unsigned long long overflow = hi(u);
		
		// Probably stall here
		__syncthreads();
		u = (BASE)u + ( TMP[base] * d) * (unsigned long long)n;
		__syncthreads();
		
		TMP[IDX] = u;
		overflow += hi(u);
		// Probably stall here
		__syncthreads();
		u = overflow + TMP[up] * ( mod(IDX , r) != (r - 1) );
		__syncthreads();
	}
	
	index++;
	base  = mod(IDX , r) > 0;

	// Propogate all carries (Rarely happens)
	do{
		__syncthreads();
		cond[index] = false;
		TMP[up] = hi( u );
		__syncthreads();
		u = (BASE)u + TMP[IDX] * base;
		if( hi(u) ) cond[index] = true;
		__syncthreads();
		// Stall
	}while( cond[index] );
	
	__syncthreads();

	for( int i = r - 1 ; i >= 0 ; i-- ){
		cond[index * ((IDX % r) == i)] = (u - n);
		__syncthreads();
		// Stall
		if( cond[index] < 0 ) {
			break;
		}

		// No stall
		if( cond[index] ) {
			u = u + (BASE)(~n) + 1 - base;
			do{
				__syncthreads();
				cond[index] = false;
				TMP[up] = hi( u );
				__syncthreads();
				u = (BASE)u + TMP[IDX] * base;
				if( hi(u) ) cond[index] = true;
				__syncthreads();
			} while( cond[index] );
			break;
		}
	}

	__syncthreads();


	R[IDX] = u;
}

// Works!
__device__ void d_mul_concurrent( P_BASE __restrict A , P_BASE __restrict B , P_BASE __restrict R , BASE n, unsigned int r , unsigned int d ) {
	
	// Used to test
	__shared__ int cond[33];

	// Result will be stored here
	unsigned long long u = 0;

	// Cache locally
	unsigned long long a = A[IDX];
	BASE b = B[IDX];

	int index = IDX / r;
	int base = index * r;
	int up = mod(IDX + 1, r) + base;

	// Perform operations
	for( unsigned int i = 0 ; i < r ; i++ ) {		
		TMP[IDX] = b;
		
		// Stall here
		u = TMP[base + i] * a + u;	// Do the computation
		TMP[IDX] = u; // save the carry
		unsigned long long overflow = hi( u );
		
		// Probably stall here
		u = (BASE)u + ( TMP[base] * d) * (unsigned long long)n;

		TMP[IDX] = u;
		overflow += hi(u);
		// Probably stall here
		u = overflow + TMP[up] * ( mod(IDX , r) != (r - 1) );
	}
	
	index++;
	base  = mod(IDX , r) > 0;

	// Propogate all carries (Rarely happens)
	do{
		cond[index] = false;
		TMP[up] = hi( u );
		u = (BASE)u + TMP[IDX] * base;
		if( hi(u) ) cond[index] = true;
		// Stall
	}while( cond[index] );

	for( int i = r - 1 ; i >= 0 ; i-- ){
		cond[index * ((IDX % r) == i)] = (u - n);

		// Stall
		if( cond[index] < 0 ) {
			break;
		}

		// No stall
		if( cond[index] ) {
			u = u + (BASE)(~n) + 1 - base;
			do{
				cond[index] = false;
				TMP[up] = hi( u );
				u = (BASE)u + TMP[IDX] * base;
				if( hi(u) ) cond[index] = true;
			} while( cond[index] );
			break;
		}
	}
	R[IDX] = u;
}

/*
// Works! (LIES, does not handle if results is greater than or equal to N)
__device__ void d_mul( C_P_BASE __restrict__ A , C_P_BASE __restrict__ B , P_BASE __restrict__ R , BASE n, unsigned int r , unsigned int d ) {
	
	__shared__ unsigned int cond;

	// Used to test if we have to subtract by N
	__shared__ int max,min;
	max = -1;
	min = -1;

	// Result will be stored here
	unsigned long long u = 0;
	
	// Cache locally
	unsigned long long a = A[IDX];
	unsigned int b = B[IDX];
	
	// Perform operations
	for( unsigned int i = 0 ; i < r ; i++ ) {
		
		__syncthreads();
		// Everything is setup, let's go!
		// Get b!
		if( IDX == i ) cond = b;
		__syncthreads();

		u = cond * a + u;	// Do the computation
		TMP[IDX] = hi( u ); // save the carry
		u = low( u ); // remove the carry


		if( IDX == 0 ) cond = u;
		__syncthreads();

		if( cond ) {
			u = u + low( cond * d) * (unsigned long long)n;
		}

		// Add the carries
		unsigned long long tmp = TMP[IDX];
		TMP[IDX] = low(u);
		u = hi( u );
		u += tmp;
		__syncthreads();
		u += TMP[(IDX + 1) % r] * ((IDX%r) != (r - 1));

	}
	

	// Propogate all carries
	do{
		TMP[(IDX + 1) % r] = hi( u );
		u= low( u );
		cond = false;		
		__syncthreads();		
		u += TMP[IDX] * (IDX > 0);
		if( hi(u) ) cond = true;
		__syncthreads();
	}while( cond);

	if( u < n ) atomicMax( &min , IDX );
	else if( u > n ) atomicMax( &max , IDX );
	__syncthreads();

	if( max > min ) {
		u += low(~n);
		u += (IDX == 0);
		do{ 
			TMP[(IDX + 1) % r] = hi( u );
			cond = false;			
			u= low( u );
			__syncthreads();
			u += TMP[IDX] * (IDX > 0);
			if( hi( u ) ) cond = true;
			__syncthreads();
		} while( cond );
	}

	R[IDX] = low(u);
}
// */

#define DELTA (r * (blockDim.x / r))


__global__ void with_sync( P_BASE A , P_BASE B , C_P_BASE __restrict M ,  BASE r , BASE d, unsigned int muls ,  P_BASE Results) {

	// Calculate when to stop
	P_BASE max = A + muls * r;

	// Calculate starting point for this block
	// A  : [blockDim.x] [blockDim.x] [blockDim.x] [blockDim.x]...
	// BDX : 0            1           2             3
	A = A + BDX * DELTA;
	B = B + BDX * DELTA;
	Results = Results + BDX * DELTA;

	BASE n = M[mod(IDX,r)];
	
	// Do work
	while( (A + IDX) < max ) {
		d_mul_concurrent_sync( A , B , Results , n , r , d );

		// Next entry
		A = A + gridDim.x * DELTA;
		B = B + gridDim.x * DELTA;
		Results = Results + gridDim.x * DELTA;
		__syncthreads();
	}
}

__global__ void no_sync( P_BASE A , P_BASE B , C_P_BASE __restrict M ,  BASE r , BASE d, unsigned int muls ,  P_BASE Results) {
	#define DELTA (r * (blockDim.x / r))

	// No work for me!
	
	// Calculate when to stop
	P_BASE max = A + muls * r;

	// Calculate starting point for this block
	// A  : [blockDim.x] [blockDim.x] [blockDim.x] [blockDim.x]...
	// BDX : 0            1           2             3
	A = A + BDX * DELTA;
	B = B + BDX * DELTA;
	Results = Results + BDX * DELTA;

	BASE n = M[mod(IDX,r)];
	
	// Do work
	while( (A + IDX) < max ) {
		d_mul_concurrent( A , B , Results , n , r , d );

		// Next entry
		A = A + gridDim.x * DELTA;
		B = B + gridDim.x * DELTA;
		Results = Results + gridDim.x * DELTA;
	}
}



#if 1
double PCFreq = 0.0;

// used to validate the data
void validate( C_P_BASE _A , C_P_BASE _B, C_P_BASE _M , C_P_BASE Results , unsigned int r , unsigned int d , int size , std::ostream &out ) {
	
	Integer M( _M , r );
	Integer R("1");
	R <<= (r * BITLEN );
	R %= M;
	R = R.modInverse( M );
	
	Integer T;
	Integer ZERO;

	unsigned int a = 0,b = 0;

	double time = 0;
	int err_loc=-1,errs=0;	
	for( int i = 0 ; i < size ;i++ ){

		Integer A( &_A[i * r] , r );
		Integer B( &_B[i * r] , r );
		Integer Result( &Results[i * r] , r );		

		T = A.montgomery_mul( B , M , d , r );
		
		if( T != Result ) {
			errs++;err_loc = i;
		}
	}

	if( errs ) {
		
		out << std::endl;
		out << "\ta : " << a << std::endl;
		out << "\tb : " << b << std::endl;

		Integer A( &_A[err_loc * r] , r );
		Integer B( &_B[err_loc * r] , r );
		Integer R( &Results[err_loc * r] , r );
		out << "\t" << errs << "\tIncorrect" << endl;
		out << "\t" << (size - errs) << "\tCorrect" << endl << endl;
		
		out << "\tA : " << A << endl;
		out << "\tB : " << B << endl;
		out << "\td : " << d << endl;
		out << "\tr : " << r << endl;

		T = A.montgomery_mul( B , M , d , r  );
		out << "" << endl;
		out << "\tGPU : " << R << endl;
		out << "\tCPU : " << T << endl << endl;
		out << "\tLast : " << T.getRawPointer()[0] << std::endl;
		
	}
}

#endif

extern "C" void run_kernel( unsigned int do_verify , int blocks , int threads , C_P_BASE data , unsigned int data_length , C_P_BASE N , unsigned int r , unsigned int d , unsigned int runs , unsigned muls_per_blocks , std::ostream &out) {

	// Output and Random Data
	out.precision(16);
	std::random_device rd;
	std::mt19937 gen(rd());

	// How much actual data
	unsigned int gpu_words = data_length >> 1;
	unsigned int gpu_bytes = sizeof( BASE ) * (gpu_words);
	unsigned int multiplications = gpu_words / r;
	
	// Allocate data on CPU
	C_P_BASE A = data;
	C_P_BASE B = &data[gpu_words];
	P_BASE R = (P_BASE)calloc( gpu_words , sizeof(BASE) ); // Results go here...

	if( R == 0 ) {
		std::cerr << "Could not allocate memory for results on CPU." << endl;
		return;
	}

	// Start Cuda
	cudaFree( 0 );
	
	// Timing
	cudaEvent_t start,stop;
	

	// Variables
	P_BASE dev_A;
	P_BASE dev_B;
	P_BASE dev_R;
	P_BASE dev_N;

	// Allocate GPU memory
	
	errorTestMsg( "\tMAL : A" , cudaMalloc( (void**)&dev_A  , gpu_bytes ) , std::cerr );
	if( dev_A == 0 ) {
		std::cerr << "\tError : Memory not allocated for A on GPU." << endl;
		free( R );
		return;
	}


	errorTestMsg( "\tMAL : B" , cudaMalloc( (void**)&dev_B  , gpu_bytes ) , std::cerr );
	if( dev_B == 0 ) {
		std::cerr << "\tError : Memory not allocated for B." << endl;
		free( R );
		cudaFree( dev_A );
		return;
	}

	errorTestMsg( "\tMAL : R" , cudaMalloc( (void**)&dev_R  , gpu_bytes ) , std::cerr );
	if( dev_R == 0 ) {
		std::cerr << "\tError : Memory not allocated for results on GPU." << endl;
		free( R );
		cudaFree( dev_A );
		cudaFree( dev_B );
		return;
	}

	errorTestMsg( "\tMAL : N" , cudaMalloc( (void**)&dev_N , sizeof(BASE) * r ) , std::cerr );
	if( dev_R == 0 ) {
		std::cerr << "\tError : Memory not allocated for N on GPU." << endl;
		free( R );
		cudaFree( dev_A );
		cudaFree( dev_B );
		cudaFree( dev_R );
		return;
	}

	cudaFuncSetCacheConfig( with_sync , cudaFuncCachePreferShared );
	cudaFuncSetCacheConfig( no_sync , cudaFuncCachePreferShared );


	// Copy data over
	errorTestMsg( "\tMCPY : A" ,	cudaMemcpy( dev_A , A , gpu_bytes		, cudaMemcpyHostToDevice ) , std::cerr );	
	errorTestMsg( "\tMCPY : B" ,	cudaMemcpy( dev_B , B , gpu_bytes		, cudaMemcpyHostToDevice ) , std::cerr );	
	errorTestMsg( "\tMCPY : N" ,	cudaMemcpy( dev_N , N , sizeof(BASE)*r	, cudaMemcpyHostToDevice ) , std::cerr );

	// Memset
	errorTestMsg( "\tMEMSET : R" ,	cudaMemset( dev_R , 0 , gpu_bytes ) , std::cerr );

	// Stats
	double gpu_sum = 0;
	float gpu_f = 0;

	unsigned int total = threads * sizeof(BASE);
	
	errorTestMsg( "\tCreate Start Event" , cudaEventCreate( &start ) , std::cerr );
	errorTestMsg( "\tCreate End Event" , cudaEventCreate( &stop ) , std::cerr );
	//cudaDeviceSetSharedMemConfig( cudaSharedMemConfig::cudaSharedMemBankSizeDefault );
	//cudaDeviceSetSharedMemConfig( cudaSharedMemConfig::cudaSharedMemBankSizeFourByte );
	//cudaDeviceSetSharedMemConfig( cudaSharedMemConfig::cudaSharedMemBankSizeEightByte );
	
	cudaDeviceSynchronize();

	// Run 
	for( int i = 0 ; i < runs ; i++ ) {
		gpu_f = 0.0f;
	
		cudaEventRecord( start , 0);
		if( threads > 32 && (32 % r) != 0) {
			with_sync<<<blocks,threads,total>>>( dev_A , dev_B , dev_N , r , d , multiplications, dev_R );
		}else {
			no_sync<<<blocks,threads,total>>>( dev_A , dev_B , dev_N , r , d , multiplications, dev_R );
		}
		errorTestMsg( "Record Event : Stop " , cudaEventRecord(stop , 0 ) , std::cerr );
		cudaEventSynchronize( stop );
		cudaEventElapsedTime( &gpu_f , start, stop );
		gpu_sum += gpu_f;
	
	}
	errorTestMsg( "\tDestroy Start Event" , cudaEventDestroy( start ) , std::cerr );
	errorTestMsg( "\tDestroy Stop Event" , cudaEventDestroy( stop ) , std::cerr );

	gpu_sum /= 1E3;


	errorTestMsg( "\tSYNC : Before" , cudaDeviceSynchronize() , std::cerr );
	errorTestMsg( "\tMCPY : R" , cudaMemcpy( R , dev_R  , gpu_bytes , cudaMemcpyDeviceToHost ) , std::cerr );
	errorTestMsg( "\tSYNC : After" , cudaDeviceSynchronize() , std::cerr );
	
	
		// Output statistics
	double gpu_seconds = gpu_sum / runs;

	out << ( gpu_seconds );
	
	out.flush();

	if( do_verify )
		validate( A , B , N  , R , r , d ,  multiplications , out );

	out.flush();


	// Shutdown Cuda
	cudaFree( dev_A );
	cudaFree( dev_B );
	cudaFree( dev_R );
	cudaFree( dev_N );

	free( R );

	cudaDeviceReset();
}
