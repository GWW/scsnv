/*
Copyright (c) 2018-2020 Gavin W. Wilson

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
#include <cstdint>
namespace gwsc {

static constexpr const char reverse_cmpl_[] = {
    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
    16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
    32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,
    48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,
    64,'T','V','G','H','E','F','C','D','I','J','M','L','K','N','O',
    'P','Q','Y','S','A','A','B','W','X','R','Z',91,92,93,94,95,
    64,'t','v','g','h','e','f','c','d','i','j','m','l','k','n','o',
    'p','q','y','s','a','a','b','w','x','r','z',123,124,125,126,127
};

class Sequence{
    public:

        const std::string & data() const { 
            return data_;
        }

        std::string & data() { 
            return data_;
        }

        void clear() {
            data_.clear();
        }

        void push_back(char c){
            data_.push_back(c);
        }

        size_t size() const {
            return data_.size();
        }

        std::string::const_iterator begin() const {
            return data_.begin();
        }

        std::string::const_iterator end() const {
            return data_.end();
        }

        template< class InputIt >
        Sequence& append( InputIt first, InputIt last ) {
            data_.append(first, last);
            return *this;
        }

        std::string::value_type operator[](size_t i) const {
            return data_[i];
        }

        std::string::value_type & operator[](size_t i) {
            return data_[i];
        }

        std::string::value_type at(size_t i) const {
            return data_[i];
        }

        std::string::value_type & at(size_t i) {
            return data_[i];
        }

        std::string::value_type rat(size_t i) const {
            return reverse_cmpl_[static_cast<unsigned int>(data_[i])];
        }

        std::string & str() {
            return data_;
        }

        const std::string & str() const {
            return data_;
        }

        void reverse_cmpl() {
	    for(size_t i = 0, j = (size() - 1); i < (size() / 2); i++, j--){
                auto tmp = rat(j);
                at(j) = rat(i);
                at(i) = tmp;
            }
            if((size() & 1) == 1){
                size_t mid = size() / 2;
                at(mid) = rat(mid);
            }
        }

    private:
        std::string data_;
};

struct ADNA4 {
    static const std::size_t size_ = 2;
    static const char * alphabet_str_;
    static const char ltable_[256];
    static const char rtable_[4];
};


struct ADNA5 {
    static const std::size_t size_ = 4;
    static const char * alphabet_str_;
    static const char   ltable_[256];
    static const char   rtable_[5];
};

template<typename T, typename R>
inline bool seq2int(const std::string & seq, R & code){
    code = 0;
    assert(seq.size() <= sizeof(R) * 8 / T::size_);
    for(size_t i = 0; i < seq.size(); i++){
        char c = T::ltable_[static_cast<size_t>(seq[i])];
        if(c < 0) return false;
        code = (code << T::size_) | c;
    }
    return code;
}

template<typename T, typename R>
inline uint32_t seq2int(const std::string & seq){
    uint32_t code = 0;
    assert(seq.size() <= sizeof(R) * 8 / T::size_);
    for(size_t i = 0; i < seq.size(); i++){
        char c = T::ltable_[static_cast<size_t>(seq[i])];
        code = (code << T::size_) | c;
    }
    return code;
}

template<typename T, typename R>
inline std::string int2seq(R val, unsigned int N){
    std::string seq("", N);
    //cout << "val = " << val << "\n";
    R mask = (1 << T::size_) - 1;
    for(size_t i = 0; i < N; i++){
        //cout << "  i = " << i << " val = " << val << " " << (val & 0x3) << "\n";
        seq[N - i - 1] = T::alphabet_str_[val & mask];
        val >>= T::size_;
    }
    return seq;
}

template <typename U>
inline static U getmask(std::size_t index, std::size_t size){
    return static_cast<U>(((1UL << (size)) - 1UL) << (index * size));
}

}
