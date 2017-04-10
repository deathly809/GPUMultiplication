#include "Scanner.h"

Scanner::Scanner(std::istream& in) {
input = &in;
}

bool Scanner::hasNext() {
	return true;
}

Integer Scanner::nextInteger() {
	return Integer();
}

unsigned int Scanner::nextUINT32() {
	return 5;
}

Scanner::~Scanner(void)
{
}
