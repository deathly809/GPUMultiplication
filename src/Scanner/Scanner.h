#pragma once

#include <istream>

#include "../Integer/Integer.h";

class Scanner
{
private:
	char buffer[1024];
	std::istream* input;
	Scanner(){};
public:
	Scanner(std::istream& input);
	~Scanner(void);

	bool hasNext();
	Integer nextInteger();
	
	unsigned int nextUINT32();

};

