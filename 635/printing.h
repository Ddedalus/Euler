#ifndef __PRINTING_H__
#define __PRINTING_H__

#include <map>

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v)
{
    os << "[";
    for (const T &p : v)
        os << p << ", ";
    os << "]";
    return os;
}

template <class K, class V>
std::ostream &operator<<(std::ostream &os, const std::unordered_map<K, V> &m)
{
    os << "{";
    for (auto p = m.begin(); p != m.end(); ++p)
        os << p->first << ":" << p->second << ", ";
    os << "}";
    return os;
};

template <class K, class V>
std::ostream &operator<<(std::ostream &os, const std::map<K, V> &m)
{
    os << "{";
    for (auto p = m.begin(); p != m.end(); ++p)
        os << p->first << ":" << p->second << ", ";
    os << "}";
    return os;
};

#endif