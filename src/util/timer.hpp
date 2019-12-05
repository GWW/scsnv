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

#include <chrono>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>

namespace gwsc {

class SimpleTimer {
    public:
        SimpleTimer()
            : t0(std::chrono::high_resolution_clock::now())
        {
        }

        void reset(){
            t0 = std::chrono::high_resolution_clock::now();
        }

        std::string elapsed() const {
            auto  t1 = std::chrono::high_resolution_clock::now();
            auto delta = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0);
            auto hours = std::chrono::duration_cast<std::chrono::hours>(delta);
            std::stringstream ss;
            ss << std::setfill('0') << std::setw(2) << hours.count() << ":";
            delta -= hours;
            auto minutes = std::chrono::duration_cast<std::chrono::minutes>(delta);
            ss << std::setw(2) << minutes.count() << ":";
            delta -= minutes;
            auto seconds = std::chrono::duration_cast<std::chrono::seconds>(delta);
            ss << std::setw(2) << seconds.count();
            delta -= seconds;
            return ss.str();
        }

        size_t seconds() const {
            auto  t1 = std::chrono::high_resolution_clock::now();
            auto delta = std::chrono::duration_cast<std::chrono::seconds>(t1-t0);
            //auto seconds = std::chrono::duration_cast<std::chrono::seconds>(delta);
            return delta.count();
        }

    protected:
        std::chrono::high_resolution_clock::time_point t0;

};

class ScopedTimer : public SimpleTimer {
    public:
        ScopedTimer(std::function<void(std::string)> callback)
            : SimpleTimer(), cb(callback)
        {
        }

        ~ScopedTimer(void)
        {
            cb(elapsed());
        }

    protected:
        std::function<void(std::string)> cb;
};

}
