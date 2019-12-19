/*
 * range.h
 *
 *  Created on: Dec 1, 2019
 *      Author: Nikhil
 */

#ifndef INCLUDE_MDSTRESSLAB___RANGE_H_
#define INCLUDE_MDSTRESSLAB___RANGE_H_

template <class T> class range {
private:
    class iter {
    private:
        T at;
    public:
        iter(T at) : at(at) {}
        bool operator!=(iter const& other) const { return at != other.at; }
        T const& operator*() const { return at; }
        iter& operator++() { ++at; return *this; }
    };

    T begin_val;
    T end_val;
public:
    range(T begin_val, T end_val) :
        begin_val(begin_val), end_val(end_val) { }
    iter begin() { return iter(begin_val); }
    iter end() { return iter(end_val); }
};




#endif /* INCLUDE_MDSTRESSLAB___RANGE_H_ */
