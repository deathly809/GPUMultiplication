#include "../Integer/Integer.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <cctype>
#include <algorithm> 
#include <functional> 
#include <locale>

// Cuda
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>	

using namespace std;

extern "C" void run_kernel( unsigned int verify ,  int blocks , int threads , C_P_BASE data , unsigned int data_len , C_P_BASE N , unsigned int r , unsigned int d , unsigned int runs , unsigned int muls_per_blocks , std::ostream &out);
void printAndExit( string msg , int return_code );
void cuda( unsigned int verify , unsigned int min_threads , unsigned int max_threads , unsigned int block_loop_start , unsigned int block_loop_end , unsigned int runs , std::istream& in , std::ostream& out );

// d = ( ( ( ( (N + 1) / 2 ) ^ BITLEN) mod N) * (2 ^BITLEN) ) / N
unsigned int compute_d( Integer& N ) {
	Integer R( N ),ONE("1");

	R += ONE;
	R >>= 1;
	R = R.exp( BITLEN , N );
	R <<= BITLEN;
	R -= ONE;
	R /= N;
	return low(R[0]);
}

int total_read = 0;

unsigned int readUINT32( std::istream& input ) {
	total_read += sizeof( unsigned int );
	union INT_CHAR {
		char C[sizeof(unsigned int)];
		unsigned int I;
	} T;
	input.read( T.C , sizeof(unsigned int) );

	
	if( input.eof() ) {
		std::cout << "eof" << std::endl;
		return 0 ;
	}
	if( input.fail() ) {
		std::cout << "failbit" << std::endl;
		return 0;
	}
		if( input.bad() ) {
		std::cout << "badbit" << std::endl;
		return 0;
	}

	return T.I;
}

unsigned int readInteger( std::istream& input , P_BASE array ) {
	unsigned int len = readUINT32( input );
	for( int i = 0 ; i < len ; i++ ) {
		BASE tmp = readUINT32( input );
		array[i] = tmp;
	}
	array[len] = 0;
	return len;
} 

void cuda( unsigned int verify , unsigned int min_threads , unsigned int max_threads , unsigned int block_loop_start , unsigned int block_loop_end , unsigned int runs , std::istream& in , std::ostream& out ) {
	// Support up to 4096-bit integers
	BASE buffer[512];

	std::string line;
	unsigned int r,count,data_len;
	
	// Get the prime number
	r = readInteger( in , buffer );
	Integer N( buffer , r );

	// Get number of integers
	count = readUINT32( in );
	
	// Allocate date
	data_len = count * r;
	P_BASE data = (P_BASE)calloc( data_len ,sizeof(BASE) );
	
	if( data == 0 ) {
		std::cerr << "Main.cpp(cuda): Could not allocated memory for data." << std::endl;
		return;
	}
	
	// Read in the numbers
	out << "Reading " << count << " integers from input" << std::endl << std::endl;
	int per = 0;
	for( int i = 0 ;  i < count ; i++ ) {
		if( (100 * i) / count >= per ) {
			out << per << "%\t";
			per += 10;
		}
		unsigned int len = readInteger( in , buffer );
		memcpy( &data[i * r] , buffer , len * sizeof( BASE ) );
	}

	unsigned int muls = count / 2;	
	unsigned int threads = max(1,(min_threads/r) ) * r;
	
	// Compute stuff needed for montgomery
	unsigned int d = compute_d( N );

	out << std::endl << std::endl << "Data" << std::endl << std::endl;

	out << "\tPrime : " << N << std::endl;
	out << "\tComputational constant r : " << r << std::endl;
	out << "\tComputational constant d : " << d << std::endl;
	out << "\tNumber of runs : " << runs << std::endl;

	unsigned int thrd = threads / r;
	while( threads <= max_threads ) {
		unsigned int start = block_loop_start;

		out << threads << "\t";

		while( start <= block_loop_end ) {
			// Number of blocks
			unsigned int blocks = (muls + thrd - 1) / thrd;
			blocks = (blocks + start - 1) / start;
			run_kernel( verify , blocks , threads , data , data_len , N.getRawPointer(), r , d  , runs , start , out );

			out << "\t";

			start++;
		}
		out << std::endl;
		thrd++;
		threads +=r;
	}

	// Free memory
	std::free( data );
}

bool is_number( const std::string& s ) {
	if( s.length() == 0 ) return false;

	std::string::const_iterator iter = s.begin();
	while( iter != s.end() ) {
		if( !std::isdigit( *iter ) ) {
			return false;
		}
		iter++;
	}
	return true;
}

void printAndExit( string msg , int return_code ) {
	std::cerr << msg << std::endl;
	//system( "PAUSE" );
	exit( return_code );
}

unsigned int tryToGet( std::map<string,unsigned int> k_v , string key  ) {
	std::map<string,unsigned int>::iterator iter = k_v.find( key );
	if( iter == k_v.end() ) {
		printAndExit( string( "Could not find flag " ) + key , EXIT_FAILURE );
	}
	return k_v[key];
}

#define testOrSet(k_v,key,value) { \
	std::map<string,unsigned int>::iterator iter = k_v.find( key ); \
	if( iter == k_v.end() ) { \
		 k_v.insert( std::pair<string,unsigned int>( key , value ) ); \
	}else { \
		value = k_v[key]; \
	} \
}	\

int main( int argc , char** args ) {

	// Holds configuration
	std::map<string,unsigned int> flags;

	// Configuration variables
	std::string		config	=	"../data/runtime.conf";
	std::string		input	=	"";
	std::string		output	=	"";
	std::istream*	in		=	&std::cin;
	std::ostream*	out		=	&std::cout;

	// Runtime parameters
	unsigned int block_loop_start	=	0;
	unsigned int runs				=	0;
	unsigned int min_threads		=	0;
	unsigned int max_threads		=	0;
	unsigned int verify				=	0;
	unsigned int device				=	0;


	// Read command line arguments
	for( int i = 1 ; i < argc ; i++ ) {
		if( args[i][0] != '-' ){
			printAndExit( string( "Invalid flag : " ) + string( args[i] ) , EXIT_FAILURE );
		}
		string s( args[i] );
		int pos = s.find_first_of( '=' );
		string key = s.substr( 1 , pos - 1 );
		
		string value = s.substr( pos + 1 );
		if( key == "config" ) {
			config = value;
		}else if( key == "input") {
			input = value;
		}else if( key == "out") {
			output = value;
		}else{
			if( is_number( value ) ) {
				std::map<string,unsigned int>::iterator it = flags.find( key );
				if( it == flags.end() ) {
					flags.insert( std::pair<string,unsigned int>( key , strtol( value.c_str() , 0 , 10 ) ) );
				}
			}else {
				std::cerr << "Expected an integer, found " << value << std::endl;
				return EXIT_FAILURE;
			}
		}
	}

	// Read in runtime 	
	in = new ifstream( config , std::ios::binary );
	if( !in->good() ) {
		printAndExit( std::string( "Could not open configuration file : " ) + config , EXIT_FAILURE );
	}

	std::string line;
	std::getline( *in , line );
	unsigned int line_number = 0;
	do {
		line = line.substr( 0 , line.size() - 1 );

		if( line.length() > 0 )  {
			if( line.substr( 0 , 1 ) != "#" ) {
				int pos = line.find( "=" );
				if( pos != string::npos ) {
					std::string key = line.substr( 0 , pos );
					std::string value = line.substr( pos + 1 );
					if( key == "input" ) {
						if( input.length() == 0) {
							input = value;
						}
					}else if( key == "out") {
						if( output.length() == 0 ) {
							output = value;
						}
					}else {
						if( is_number( value ) ) {
							std::map<string,unsigned int>::iterator it = flags.find( key );
							if( it == flags.end() ) {
								flags.insert( std::pair<string,int>( key , strtol( value.c_str() , 0 , 10 ) ) );
							}
						}else {
							printAndExit( std::string( "Expected an integer got : " ) + value + std::string( " on line " ) + std::to_string( (long long)line_number) , EXIT_FAILURE );
						}
					}
				} else {
					printAndExit( std::string("Error on line number " ) + std::to_string( (long long)line_number ) + std::string( "\n" ) + line , EXIT_FAILURE );
				}
			}
		}
		std::getline( *in , line );
		line_number++;
	}while( !in->eof() );

	// multiplications per block
	block_loop_start = 1;
	testOrSet( flags , "block_loop_start" , block_loop_start );
	
	int block_loop_end = block_loop_start;
	testOrSet( flags , "block_loop_end" , block_loop_end );
	
	// Threads
	min_threads = 1;
	testOrSet( flags , "min_threads" , min_threads );
	
	max_threads = min_threads;
	testOrSet( flags , "max_threads" , max_threads );

	// runs
	runs = 1;
	testOrSet( flags , "runs" , runs );
	
	// Verify
	verify = 0;
	testOrSet( flags , "verify" , verify );

	// Input
	if( input.length() != 0 ) {
		in = new ifstream( input , std::ios::in | std::ios::binary );
	}else{
		*out << "Input set to cin" << std::endl;
		in = &cin;
	}
	
	// Output
	if( output.length() != 0 ) {
		out = new ofstream( output );
	}else{
		out = &cout;
	}

	// Select Device
	device = 0;
	testOrSet( flags , "device" , device );

	cudaSetDevice( device );

	cudaDeviceProp devProp;
	cudaGetDeviceProperties( &devProp , device );
	*out << "Using Device : " << devProp.name << std::endl;

	// Run
	cuda( verify , min_threads , max_threads , block_loop_start , block_loop_end, runs , *in , *out );

	if( in != &cin ) {
		delete( in );
	}

	if( out != &cout ) {
		delete( out );
	}

	//system("PAUSE");

	return EXIT_SUCCESS;
}
