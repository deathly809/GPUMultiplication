#pragma once

/*
*	Author : Jeffrey A. Robinson
*  Date   : August 25, 2013
*  "Arbitrary" Positive Integer library.  Only works with positive integers and expected to be used with number theory applications.
*/

// need if we do DLL
#ifdef WINDOWS
#define INTEGER_API __declspec( dllexport )
#else
	#define INTEGER_API
#endif

#include <ostream>
#include <stdlib.h>
#include <sstream>

#include <string>

#define BASE unsigned int
#define P_BASE BASE*
#define C_P_BASE const P_BASE

#define BITLEN (sizeof(BASE) * 8)
#define MASK (unsigned int)( ( ( (1 << (BITLEN - 1) ) - 1 ) << 1 ) + 1 )

#define low(x) ((x) & MASK)
#define hi(x) (((x)>>BITLEN) & MASK)
#define split(x,h,l) { h = (BASE)hi(x) ; l = (BASE)low(x); }
#define combine(x,h,l) { x = unsigned long long(low(h) << BITLEN) | unsigned long long(low(l)); }
#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

class INTEGER_API Integer {
private:
	BASE *data;
	BASE len;
	char sng;

	void compress();	
	void read( const char* , unsigned int );

public:
	Integer();
	Integer( unsigned int size );
	Integer( const char* string );
	Integer( char* string , unsigned int len );
	Integer( const Integer& other );
	Integer( C_P_BASE data , int len );
	
	~Integer();

	// Stuff
	unsigned int getWordCount();
	C_P_BASE getRawPointer() { return this->data; }
	void setLength( unsigned int len );
	bool isZero();

	// Math
	Integer& operator+=(const Integer& rhs);
	Integer& operator-=(const Integer& rhs);
	Integer& operator*=(const Integer& rhs);
	Integer& operator/=(const Integer& rhs);
	Integer& operator%=(const Integer& rhs);
	Integer& operator>>=(unsigned int);
	Integer& operator<<=(unsigned int);
	Integer& operator++(int);
	Integer& operator--(int);

	// Access
	BASE operator[] ( const int x ) const;
	Integer& operator=(const Integer& other);

	// Comparison
	bool operator==(const Integer& other) const;
	bool operator!=(const Integer& other) const;
	bool operator<(const Integer& other) const;
	bool operator<=(const Integer& other) const;
	bool operator>(const Integer& other) const;
	bool operator>=(const Integer& other) const;

	// Special Math!
	Integer montgomery_mul( const Integer& other , const Integer& N , unsigned int d , unsigned int r );
	Integer modInverse( const Integer& N );
	void div( const Integer& rhs , Integer& quotent, Integer& remainder );
	Integer exp( unsigned int , const Integer& N);
	Integer exp( const Integer &pow , const Integer& N );
	Integer gcd( const Integer & other );

	friend std::ostream& operator<<(std::ostream& os,const Integer& op);
	friend Integer operator%( const Integer& lhs , const Integer& rhs );
	friend Integer operator+( const Integer& lhs , const Integer& rhs );
	friend Integer operator-( const Integer& lhs , const Integer& rhs );
	friend Integer operator/( const Integer& lhs , const Integer& rhs );
	friend Integer operator*( const Integer& lhs , const Integer& rhs );
	friend Integer operator<<( const Integer& h , unsigned int);
	friend Integer operator>>( const Integer& h ,unsigned int);
};


inline Integer operator+( const Integer& lhs , const Integer& rhs ) {
	Integer res = lhs;
	res += rhs;
	return res;
}
inline Integer operator-( const Integer& lhs , const Integer& rhs ) {
	if( &lhs == &rhs ) return Integer();
	Integer res = lhs;
	res -= rhs;
	return res;
}

inline Integer operator*( const Integer& lhs , const Integer& rhs ) {
	Integer res = lhs;
	res *= rhs;
	return res;
}
inline Integer operator/( const Integer& lhs , const Integer& rhs ) {
	Integer res = lhs;
	res /= rhs;
	return res;
}

inline Integer operator%( const Integer& lhs , const Integer& rhs) {
	Integer res = lhs;
	res %= rhs;
	return res;
}

inline Integer operator<<(const Integer& lhs , unsigned int shift ) {
	Integer res = lhs;
	res <<= shift;
	return res;
}

inline Integer operator>>(const Integer& lhs , unsigned int shift ) {
	Integer res = lhs;
	res >>= shift;
	return res;
}

inline std::ostream& operator<<(std::ostream& out , const Integer& in ){

	Integer ZERO;

	if( in == ZERO ) {
		out << "0";
		return out;
	}

	Integer DIV( "10" );
	unsigned long long l = 0;

	std::string output = "";

	std::stringstream str;

	Integer t1(in);
	Integer t2;
	int pos = in.len * 20 - 1;

	while( ( t1 > ZERO ) && pos >= 0) {
		t1.div( DIV , t1 , t2 );
		l = t2.data[1];
		l <<= BITLEN;
		l |= t2.data[0];

		output = std::to_string( l ) + output;

		pos--;
	}

	//str >> output;

	if( pos < 0 ) {
		out << "There was an Error printing the number." << std::endl;
		return out;
	}
	
	out << output;
	return out;
}
