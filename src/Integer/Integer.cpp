/*
*	Author : Jeffrey A. Robinson
*  Date   : August 25, 2013
*  "Arbitrary" Positive Integer library.  Only works with positive integers and expected to 
*	be used with number theory applications.
*/

#include "../Integer/Integer.h"

#include <cstring>
#include <string>
#include <iostream>

using namespace std;


// Private


void Integer::compress( ) {
	int i = this->len - 1;
	while( this->data[i] == 0 ) {
		i--;
		if( i < 0 ) break;
	}
	setLength( max( 1 , i + 1 ) );
}

void Integer::read( const char* string , unsigned int len ) {
	Integer res(1);
	Integer other(1);
	Integer TEN(1);
	TEN.data[0] = 10;

	unsigned int i = 0;
	while( i < len ) {
		other.data[0] = string[i] - '0';
		res *= TEN;
		res += other;
		i++;
	}
	setLength( res.len );
	memcpy( data , res.data , res.len * sizeof( BASE ) );
}

unsigned int Integer::getWordCount() { 
	compress();
	return this->len; 
}

// Constructor
Integer::Integer() : len(0), data(0) { 
	setLength( 1 );
}
Integer::Integer( unsigned int size ): len(0), data(0) { 
	setLength( size );
}
Integer::Integer( const char* string ): len(0), data(0) {
	setLength( 1 );
	read( string , strlen( string ) );
	compress();
}
Integer::Integer( char* string , unsigned int len  ): len(0) , data(0){
	setLength( 1 );
	read( string , len );
}
Integer::Integer( const Integer& other ) : len(0), data(0){
	if( &other == this ) return;

	setLength( other.len );

	if( other.data && data ) {
		memcpy( data , other.data , sizeof( BASE) * (other.len) );
	}
	compress();
}
Integer::Integer( C_P_BASE data , int len ) : len(0), data(0) {
	setLength( len );
	if( data && this->data ) {
		memcpy( this->data , data , len * sizeof( BASE ) );
	}
	compress();
}
// Deconstructor
Integer::~Integer() { 
	free( data );
}

bool Integer::isZero() {
	for( int i = 0 ; i < this->len ; i++ ) {
		if( this->data[i] ) return false;
	}
	return true;
}


void Integer::setLength( unsigned int new_len ) {
	if( this->len != new_len ) {
		P_BASE d = (P_BASE)realloc( this->data , sizeof( BASE) * (new_len + 1) );
		if( d != NULL ) {
			this->data = d;
		}
	}
	for( unsigned int i = this->len ; i <= new_len ; i++ ) {
		this->data[i] = 0;
	}	
	this->len = new_len;
}

// Math

Integer& Integer::operator+=(const Integer& rhs){
	Integer ZERO;
	Integer lhs = *this;
	
	if( (rhs == ZERO) ) {
		return *this;
	}if( lhs == ZERO ) {
		*this = rhs;
		return *this;
	}

	if( lhs.len < rhs.len ) {
		lhs.setLength( rhs.len );
	}

	unsigned int size = 0;

	P_BASE res = (P_BASE)calloc( lhs.len + 1 , sizeof(BASE) );
	if( res == 0 ) {
		std::cerr << "Integer.cpp(+=): Could not allocated memory for res." << std::endl;
		return *this;
	}

	// Add lower
	unsigned int carry = 0;
	for( size = 0 ; size < rhs.len; size++ ) {
		unsigned long long tmp = lhs.data[size];
		tmp += rhs.data[size];
		tmp += carry;
		split( tmp , carry , res[size] );
	}

	// Carry higher
	for( ; size < lhs.len ; size++ ) {
		unsigned long long tmp = lhs.data[size];
		tmp += carry;
		split( tmp , carry , res[size] );
	}res[size] = carry;

	setLength( lhs.len + 1 );
	memcpy( data , res , ( len ) * sizeof( BASE ) ) ;
	free( res );

	compress();
	return *this;
}

Integer& Integer::operator-=(const Integer& rhs ) {
	if( this == &rhs ) {
		setLength(1);
		data[0] = 0;
		return *this;
	}
	
	Integer ZERO;

	// Test special cases
	if( rhs == ZERO ) {
		*this = rhs;
		return *this;
	}
	
	// BUG: This is not correct correct
	if( *this < rhs ) {
		*this = (rhs - *this);
		return *this;
	}
	
	// Create Return value
	P_BASE res = (P_BASE)calloc( len + 1, sizeof( BASE ) );
	if( res == 0 ) {
		std::cerr << "Integer.cpp(-=): Could not allocated memory for res." << std::endl;
		return *this;
	}
	
	memcpy( res , rhs.data , sizeof( BASE ) * rhs.len );

	// Perform subtraction
	unsigned int carry = 1;
	for( unsigned int i = 0 ; i < len ; i++ ) {
		unsigned long long tmp = low(~(res[i]));
		tmp += carry;
		tmp += data[i];
		split( tmp , carry , res[i] );
	}

	memcpy( data , res , sizeof( BASE ) * len );
	data[len] = 0;

	free( res );

	// Compress the results
	compress();
	return *this;
}

Integer& Integer::operator*=(const Integer& rhs) {
	std::size_t size = len + rhs.len + 1;

	P_BASE res = (P_BASE)calloc( size , sizeof(BASE) );	
	if( res == 0 ) {
		std::cerr << "Integer.cpp(*=): Could not allocated memory for res." << std::endl;
		return *this;
	}
	unsigned int carry = 0;
	for( unsigned long long i = 0 ; i < this->len ; i++ ) {
		carry = 0;
		int k = i;
		for( unsigned long long j = 0 ; j < rhs.len ; j++ ) {
			unsigned long long prod = this->data[i];
			prod *= rhs.data[j];
			prod += res[k];
			prod += carry;
			split( prod , carry , res[k] );
			k++;
		}
		res[k] = carry;
	}
	setLength( size );
	memcpy( data , res , size * sizeof( BASE ) );
	free( res );

	compress();
	return *this;
}
Integer& Integer::operator%=(const Integer& rhs ) {
	Integer q,r;
	this->div( rhs , q , r );
	setLength( r.len );
	memcpy( data , r.data , r.len * sizeof( BASE ) );
	compress();
	return *this;
}

Integer& Integer::operator/=(const Integer& rhs) {
	Integer q,r;
	this->div( rhs , q , r );
	setLength( q.len );
	memcpy( data , q.data , q.len * sizeof( BASE ) );
	compress();
	return *this;
}

Integer& Integer::operator++(int) {
	int carry = 1;
	for( int i = 0 ; (i < len) && carry; i++ ) {
		unsigned long long t = carry;
		t += data[i];
		split( t , carry , data[i] );
	}

	if( carry ) {
		setLength( len + 1 );
		data[len - 1] = carry;
	}
	compress();
	return *this;
}
Integer& Integer::operator--(int) {
	Integer ONE("1");
	Integer res = (*this) - ONE;
	memcpy( data , res.data , len * sizeof( BASE) ) ;
	compress();
	return *this;
}

// Overloaded Access/Assignment
BASE Integer::operator[](int pos) const {
	return this->data[pos % this->len];
}
Integer& Integer::operator=( const Integer& other) {
	int t = other.len;
	setLength( t );
	memcpy( data , other.data , t * sizeof( BASE ) );
	return *this;
}


// Overloaded Shifting
Integer& Integer::operator<<=( unsigned int shift ) {

	// Special Cases
	if( shift == 0 ) return *this;
	if( this->data[len-1] == 0 ) return *this;

	// Create result
	setLength( this->len + (shift + BITLEN - 1 ) / BITLEN  );

	// local variables
	int i = this->len;
	int d = shift / BITLEN;
	int s = shift % BITLEN;

	// Shift bits
	if( s > 0 ) {
		for( ; i ; i-- ) {
			unsigned long long tmp = this->data[i];
			tmp <<= s;
			tmp = tmp | (this->data[i - 1] >> (BITLEN - s));
			this->data[i] = low(tmp);
		}this->data[0] = this->data[0] << s;
	}
	// Shift blocks
	if( d > 0 ) {
		for( int i = this->len ; i >= d ; i-- ) {
			this->data[i] = this->data[i - d];
			this->data[i - d] = 0;
		}
	}
	compress();
	return *this;
}
Integer& Integer::operator>>=( unsigned int shift ) {

	if( shift == 0 ) return *this;

	unsigned int i;
	int d = shift / BITLEN;
	int s = shift % BITLEN;

	// Shift bits
	unsigned long long x,y;
	for( i = 0 ; i < this->len ; i++ ) {
		x = this->data[i + 1];
		y = this->data[i];
		x = (x << BITLEN) + y;
		this->data[i] = low(x >> s);
	}

	// Shift blocks
	for( i = 0 ; i < (this->len - d ); i++ ) {
		this->data[i] = this->data[i + d];
	}
	for( ; i < this->len ; i++ ) {
		this->data[i] = 0;
	}
	compress();
	return *this;
}

// Overloaded Logical Operators
bool Integer::operator==(const Integer& rhs ) const {
	if( this == &rhs ) return true;

	int comp = len - 1;
	if( len < rhs.len ) {
		comp = rhs.len - 1;
		while( comp >= len ) {
			if( rhs.data[comp] > 0 ) {
				return false;
			}
			comp--;
		}
	}else if( len > rhs.len ){
		while( comp >= rhs.len ) {
			if( data[comp] > 0 ) {
				return false;
			}
			comp--;
		}
	}
	for( int i = comp ; i >= 0 ;i-- ) {
		if( this->data[i] != rhs.data[i] ) {
			return false;
		}
	}
	return true;
}
bool Integer::operator!=(const Integer& rhs ) const  {
	return !((*this)==rhs);
}
bool Integer::operator<(const Integer& rhs ) const {
	int comp = this->len - 1;
	// rhs has the possiblility of being larger (false)
	if( len < rhs.len ) {
		comp = rhs.len - 1;
		while( comp >= len ) {
			if( rhs.data[comp] > 0 ) {
				return true;
			}
			comp--;
		}
	}else if( len > rhs.len ){ // Check high order bits!
		while( comp >= rhs.len ) {
			if( data[comp] > 0 ) {
				return false;
			}
			comp--;
		}
	}
	// Check where they match in size
	while( comp > 0) {
		if( data[comp] != rhs.data[comp] ) {
			return data[comp] < rhs.data[comp];
		}
		comp--;
	}
	return data[comp] < rhs.data[comp];
}
bool Integer::operator<=(const Integer& rhs ) const {
	return !((*this)>rhs);
}
bool Integer::operator>(const Integer& rhs ) const {
	int comp = len - 1;
	// rhs has the possiblility of being larger (false)
	if( len < rhs.len ) {
		comp = rhs.len - 1;
		while( comp >= len ) {
			if( rhs.data[comp] > 0 ) {
				return false;
			}
			comp--;
		}
	}else if(len > rhs.len ) { // Check high order bits!
		comp = len - 1;
		while( comp >= rhs.len ) {
			if( data[comp] > 0 ) {
				return true;
			}
			comp--;
		}
	}
	// Check where they match in size
	while( comp ) {
		if( data[comp] > rhs.data[comp] ) {
			return true;
		}
		comp--;
	}
	return data[comp] > rhs.data[comp];
}
bool Integer::operator>=(const Integer& rhs ) const {
	return !((*this)<rhs);
}




// Other

Integer Integer::gcd( const Integer &other ) {
	Integer a = *this;
	Integer b = other;
	
	while ( !a.isZero() ) {
		Integer c = a;
		a = b%a;
		b = c;
	}
	return b;
}

Integer Integer::modInverse( const Integer& N ) {
	Integer ONE,ZERO;
	ONE.data[0] = 1;

	Integer _new("1");
	Integer _old;
	Integer q( N );
	Integer a( *this );

	Integer r;
	Integer h;
	

	int pos = 0;
	while( a > ZERO ) {
		q.div( a , q , r );
		h = q * _new;
		h = _old + h;

		_old = _new;
		_new = h;
		q = a;
		a = r;

		pos = !pos;
	}

	if( q != ONE ) {
		return ZERO;
	}
	if( pos ) {
		return _old;
	}else {
		return N  - _old;
	}
}
void Integer::div( const Integer& rhs , Integer& quotiont , Integer& remainder ) {
	Integer lhs = *this;

	Integer ONE;
	ONE.data[0] = 1;

	if( lhs < rhs ) {
		quotiont.setLength( 1 );
		quotiont.data[0] = 0;
		remainder = lhs;
		return;
	}

	if( lhs == rhs || &lhs == &rhs ) {
		quotiont.setLength( 1 );
		quotiont.data[0] = 1;

		remainder.setLength( 1 );
		remainder.data[0] = 0;
		return;
	}

	// Do computation
	Integer r_q;
	Integer r_r;

	int i = this->len * BITLEN;
	while( i > 0) {
		i--;
		unsigned int shift = (i) % BITLEN;

		r_q <<= 1;
		r_r <<= 1;

		r_r.data[0] = r_r.data[0] | ( (data[i / BITLEN] >> shift) & 1);

		if( r_r >= rhs ) {
			r_r -=rhs;
			r_q++;
		}		
	}

	if( r_r >= rhs ) {
		r_r -=rhs;
		r_q++;
	}

	r_q.compress();
	r_r.compress();

	quotiont.setLength( r_q.len );
	memcpy( quotiont.data , r_q.data , r_q.len * sizeof( BASE ) );

	remainder.setLength( r_r.len );
	memcpy( remainder.data , r_r.data , r_r.len * sizeof( BASE ) );
}

Integer Integer::exp( unsigned int x, const Integer& N ) {
	Integer base( *this );
	Integer result( "1" );
	while( x > 0 ) {
		if( x % 2 ) {
			result *= base;
			result %= N;
		}
		x >>= 1;
		base *= base;
		base %= N;
	}
	result.compress();
	return result;
}

Integer Integer::exp( const Integer& p, const Integer& modulus ) {
	Integer ZERO;
	Integer exponent(p);
	Integer base(*this);
	Integer result("1");


	Integer R( modulus );
	R++;
	R >>= 1;
	R = R.exp( BITLEN , modulus );
	R <<= BITLEN;
	R--;
	R /= modulus;

	BASE d = low(R[0]);

	while( exponent > ZERO ) {
		if( exponent.data[0] % 2 ) {
			//result = result.montgomery_mul( base , modulus , d , modulus.len );
			result *= base;
			result %= modulus;
		}
		exponent >>= 1;
		//base = base.montgomery_mul( base , modulus , d , modulus.len );
		base *= base;
		base %= modulus;
	}

	result.compress();
	return result;
}


Integer Integer::montgomery_mul( const Integer& other , const Integer& N  , unsigned int d , unsigned int r ) {
	
	P_BASE A = (P_BASE)calloc( r , sizeof(BASE) );
	if( A == 0 ) {
		std::cerr << "montgomery_mul : Could not allocated memory for A." << std::endl;
		return Integer();
	}

	P_BASE B = (P_BASE)calloc( r , sizeof(BASE) );
	if( B == 0 ) {
		free(A);
		std::cerr << "montgomery_mul : Could not allocated memory for B." << std::endl;
		return Integer();
	}

	unsigned long long *res = (unsigned long long *)calloc( r + 1 , sizeof( unsigned long long ) );
	if( res == 0 ) {
		free(A);
		free(B);
		std::cerr << "montgomery_mul : Could not allocated memory for res." << std::endl;
		return Integer();
	}
	unsigned long long *tmp = (unsigned long long *)calloc( r + 1 , sizeof( unsigned long long ) );
	if( tmp == 0 ) {
		free(A);
		free(B);
		free( res );
		std::cerr << "montgomery_mul : Could not allocated memory for tmp." << std::endl;
		return Integer();
	}

	memcpy( A , data , len * sizeof( BASE ) );
	memcpy( B , other.data , other.len * sizeof( BASE ) );

	for( int i = 0 ; i < r ; i++ ) {

		// Multiply
		// A * B[i]
		unsigned int c,carry = 0;
		for( int j = 0 ; j < r; j++ ) {
			tmp[j] = A[j];
			tmp[j] = tmp[j] * B[i];
			split( tmp[j] , c , tmp[j] );
			tmp[j] += carry;
			carry = c;
		}tmp[r] = carry;

		carry = 0;
		for( int j = 0 ; j <= r ; j++ ) {
			res[j] = tmp[j] + res[j] + carry;
			split( res[j], carry , res[j] );
		}

		if( res[0] ) {
			// TMP = N * [ (T[0] * d) % 2^b ]
			unsigned long long t = low(res[0] * d);
			carry = 0;
			for( unsigned long long j = 0 ; j < r ; j++ ) {
				tmp[j] = (N.data[j] * t);
				split(tmp[j],c,tmp[j]);
				tmp[j] += carry;
				carry = c;
			}
			tmp[r] = carry;

			carry = 0;
			for( int j = 0 ; j <= r ; j++ ) {
				res[j] = res[j] + tmp[j] + carry;
				split( res[j] , carry , res[j] );
			}

		}

		// Shift
		for( unsigned int j = 0 ; j < r ; j++ ) {
			res[j] = res[j + 1];
		}res[r] = 0;

	}

	unsigned int carry = 0;
	Integer result( r + 1 );
	for( unsigned int i = 0 ; i <= r ; i++ ) {
		res[i] += carry;
		carry = hi( res[i] );
		result.data[i] = low(res[i]);
	}
	result.compress();

	if( result >= N ) {
		result %= N;
	}

	free( A );
	free( B );
	free( res );
	free( tmp );

	result.compress();

	return result;
}