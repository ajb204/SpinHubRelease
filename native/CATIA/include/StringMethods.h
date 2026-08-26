/*
 * StringMethods.h
 *
 *  Created on: Jun 16, 2009
 *      Author: guillaume
 */


#ifndef STRINGMETHODS_H_
#define STRINGMETHODS_H_

#include <standard.h>

#ifndef MAX_STRING_LENGTH
#define MAX_STRING_LENGTH 1024
#endif

std::vector<std::string> split(std::string, const char*);
void Tab2Space(char*, int);
void ClearBuf(char*, int);
bool IsComment(char*);
bool IsComment(const std::string&);

#endif /* STRINGMETHODS_H_ */
