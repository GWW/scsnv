/*
Copyright (c) 2018-2019 Gavin W. Wilson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#pragma once

#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>

namespace gwsc {

// Class for quickly parsing a str::string into char* tokens
// both this class and Tokenizer replace deliminators with null characters
// the whole point of this is to avoid having to copy string values around
class Tokenizer {
    public:
        typedef std::vector<char *> tokens;

        Tokenizer(std::string & str, char delim = '\t') : delim_(delim), curr_(0), col_num_(-1) {
            str.c_str();
            str_ = &str[0];
        }

        Tokenizer(char * str, char delim = '\t') : str_(str), delim_(delim), curr_(0), col_num_(-1) 
        {
        }

	static void get(std::string & str, char delim, tokens & v){
	    Tokenizer tk(str, delim); 
	    tk.get_all(v);
	}

	static void get(char * str, char delim, tokens & v){
	    Tokenizer tk(str, delim); 
	    tk.get_all(v);
	}

        void get_all(tokens & v) {
            char *col;
	    v.clear();
            while((col = next()) != NULL){
                v.push_back(col);
            }
        }

        bool has_next() const {
            return str_[curr_];
        }

        char * next() {
            start_ = curr_;
            if(!str_[curr_]) return NULL;

            while(str_[curr_] && str_[curr_] != delim_) curr_++;
            if(str_[curr_]) {
                str_[curr_] = '\0';
                curr_++;
            }
            col_num_++;
            // Move past the null character
            return str_ + start_;
        }

        size_t start() const {
            return start_;
        }
        int col_num()  const {
            return col_num_;
        }

    private:
        char       * str_;
        char         delim_;
        size_t       curr_;
        size_t       start_;
        int          col_num_;
};

}
