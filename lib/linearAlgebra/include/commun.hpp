#ifndef COMMUN_HPP
#define COMMUN_HPP

// #ifndef ARDUINO
// typedef unsigned long long size_t;
// #include <string>
// typedef std::string String;
// #endif
#include <Arduino.h>

template <typename T = float>
class Vector;

template <typename T>
class rowMajorMatrix;

template <typename T>
class symMatrix;


namespace internal
{
    // float InvSqrt43 (float x);
    static size_t vec_count = 0;
    static size_t alloc_count = 0;
    template <typename T> T &&move(T &a) noexcept;
    template <typename T> void swap(T &a, T &b) noexcept;
    template <class Derived>
    class tmp;
    template <typename T> static const T _zero = T();
    template <typename T> static const T _one = T(1);
} // namespace internal;

namespace operators
{
    static bool VectorCheckSize = true;
    static bool MatrixCheckSize = true;
} // namespace operators


#include "vector.hpp"
namespace internal
{
    static size_t tmp_count = 0;
    template <class Derived>
    class tmp : public Derived
    {
    private:
        static Vector<tmp *> buffer;

    public:
        bool currentlyUsed = true;
        using Derived::Derived;
        template <typename... Args>
        static tmp *get(Args... shape);
        tmp *release();
        static const size_t currentlyUsedCount();
        static const size_t bufferSize() { return buffer.size(); }
        static void freeAll();
    };

} // namespace internal

#ifndef COMMUN_CPP
#include "commun.cpp"
#endif // COMMUN_CPP

#endif // COMMUN_HPP