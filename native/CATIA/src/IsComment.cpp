/*
 * IsComment.cpp
 *
 *  Created on: Jun 17, 2009
 *      Author: guillaume
 */

#include <standard.h>

#include <StringMethods.h>

bool IsComment(char *line) {
	return line[0] == '#';
}

bool IsComment(const std::string& line) {
	return line[0] == '#';
}
