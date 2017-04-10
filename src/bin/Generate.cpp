
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cctype>
#include "../Integer/Integer.h"
#include <time.h>

#include <random>

int total_written = 0;
std::mt19937_64 gen;

unsigned long long rand64() {
	return gen();
}
using namespace std;

bool isNumber( string s ) {
	if( s.length() == 0 ) return 0;
	std::string::iterator iter = s.begin();
	while( iter != s.end() ) {
		if( !std::isdigit( *iter ) ) {
			return false;
		}
		iter++;
	}
	return true;
}

Integer generate( unsigned int bits ) {

	int words = (bits + BITLEN - 1) / BITLEN;
	P_BASE data = (P_BASE)calloc( words , sizeof(BASE ) );
	if( data == 0 ) {
		std::cerr << "Generate.cpp(generate): Could not allocated memory for data." << std::endl;
		return Integer();
	}

	for( int i = 0 ; i <  words - 1 ; i+=2 ) {
		split( rand64() , data[i + 1] , data[i] );
	}

	if( words % 2 ) {
		data[words - 1] = low( rand64() );
	}

	unsigned long long _MASK = 1;
	_MASK <<= (bits) % (BITLEN + 1);
	_MASK -= 1;
	data[words-1] &= _MASK;

	return Integer( data , words );
}

bool fermat_test( Integer& P ) {
	BASE TEST = P.getRawPointer()[0];	
	if( TEST % 2 == 0 ) {
		std::cout << "Even number!" << std::endl;
		return false;
	}

	Integer P_1( P );
	P_1--;

	Integer ONE( "1" );
	
	for( int i = 0 ; i < 1000 ; i++ ) {		
		Integer a = generate( (P.getWordCount() * BITLEN) - 1 );
		a = a.exp( P_1 , P );
		if( a != ONE ) {
			return false;
		}
	}
	return true;
}

void writeUINT32( unsigned int value , std::ostream& output ) {
	total_written += sizeof( value );
	union INT_CHAR {
		char C[sizeof(unsigned int)];
		unsigned int I;
	} T;
	T.I = value;
	output.write( T.C , sizeof(unsigned int) );
	output.flush();
}

void writeInteger( Integer& a , std::ostream& out ) {
	unsigned int len = a.getWordCount();
	writeUINT32( len , out );
	C_P_BASE data = a.getRawPointer();
	for( int i = 0 ; i < len ; i++ ) {
		writeUINT32( data[i] , out );
	}
	if( len != a.getWordCount() ) {
		std::cerr << "Word Count Changed!" << std::endl;
		std::cerr << "From : " << len << " To : " << a.getWordCount() << std::endl;
	}
}


int main(int argc, char** args ) {

	Integer TWO( "2" );

	// Random number generator
	gen.seed( time( NULL ) );

	// Check to make sure we have enough parameters
	if( argc != 4 ) {
		cerr << args[0] << " <prime length> <generate count> <output>" << endl;
		system( "PAUSE" );
		return EXIT_FAILURE;
	}

	// Validate some things
	if( !isNumber( args[1] ) ) { 
		cerr << "Prime number length is not a number!" << endl;
		system( "PAUSE" );
		return EXIT_FAILURE;
	}

	if( !isNumber( args[2] ) ) { 
		cerr << "Count is not a positive number!" << endl;
		system( "PAUSE" );
		return EXIT_FAILURE;
	}


	unsigned int count = std::strtoul( args[2] , 0 , 10 );
	if( count == 0 ) {
		cerr << "Count is not a positive number!" << endl;
		system( "PAUSE" );
		return EXIT_FAILURE;
	}

	unsigned int bits = std::strtoul( args[1] , 0 , 10 );
	if( count == 0 ) {
		cerr << "Length is not a positive number!" << endl;
		system( "PAUSE" );
		return EXIT_FAILURE;
	}
	
	// Try to open the file
	std::ofstream output;
	output.open( args[3] , std::ios::binary );


	if( output.good() ) {

		// Create the prime number
		Integer N;
		

		// 1 : Generate a random number
		N = generate( bits );
		// If even...
		
		if( N.getRawPointer()[0] % 2 == 0 ) {
			N++;
		}


		// 2 : Test if number is "maybe" prime
		while( !fermat_test( N ) ) {
			// 3 : If not "maybe" prime generate another number
			N += TWO;
		}
		
		// Write N to file
		writeInteger( N , output );
		// Write the number of integers generated to a file
		writeUINT32( count , output );

		
		Integer R( N ),ONE("1"),ZERO;
		R += ONE;
		R >>= 1;
		R = R.exp( BITLEN , N );
		R <<= BITLEN;
		R -= ONE;
		R /= N;
		BASE d = low(R[0]);


		// Generate a random beginning number
		BASE r[2];
		Integer out("1");	
		for( int i = 0; i < N.getWordCount(); i+=2 ) {
			split( rand64() , r[0] , r[1] );
			Integer mul(r , 2 );
			out = out.montgomery_mul( mul , N , d , N.getWordCount() );
		}

		// Generate the numbers to output
		for( int i = 0 ; i < count ; i++ ) {
			
			if( out >= N ) {
				std::cout << "Error!" << std::endl;
			}

			// Write the current number to a file
			writeInteger( out , output );
			
			// Generate the next number
			split( rand64() , r[0] , r[1] );
			Integer mul(r , 1 );
			out = out.montgomery_mul( mul , N , d , N.getWordCount() );
			out %= N;
		}

		// Close the file
		output.flush();
		output.close();
	}else {
		// well I guess we could not create the file
		std::cerr << "Could not open the file." << std::endl;
	}

	// Output the total number of bytes written

	return EXIT_SUCCESS;
}