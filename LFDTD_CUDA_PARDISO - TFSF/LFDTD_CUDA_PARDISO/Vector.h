/********************************************************************************************************************************
* Created by Yifan Wang on 5/23/24.
*
* Currently supported operation:
*
************************************************  Level I BLAS  *****************************************************************
* = : let this = RHS vector (same type)                            update current vector
* = (overloaded): let this = RHS vector (different type)           update current vector
* + : V = this + RHS vector                                        return a new vector
* + (overloaded): V = this + const_val                             return a new vector
* += : this += const_val                                           update current vector
* += (overloaded) : this += RHS vector                             update current vector
* - : V = this - RHS vector                                        return a new vector
* - (overloaded): V = this - const_val                             return a new vector
* -= : this -= const_val                                           update current vector
* -= (overloaded): this -= RHS vector                              update current vector
* * : V = this * RHS vector                                        return a new dot product result regardless dimension type
* * (overloaded): V = this * const_val                             return a new vector
* *= : this * = const_val                                          update current vector
*********************************************************************************************************************************
*
*
************************************************  Common Operator  **************************************************************
* operator [int idx] : return boost::multi_array::reference
* operator [int idx] const : return boost::multi_array::const_reference
* operator (int idx) : return this vector value at index idx regardless dimension type - Deprecated
* operator (int idx) const : return const this vector value at index idx regardless dimension type - Deprecated
* operator == : return true if this == RHS in terms of m and n length, regardless of the type
* operator != : return true if this != RHS in terms of m and n length, regardless of the type
*********************************************************************************************************************************
*
*
* ************************************************  Member Function  ************************************************************
* constexpr void fill_v(const Type& Val) - fill in the vec with specific value
* Type M() const - return size M
* Type N() const - return size N
* template<typename A> bool check_integrity(const vec<A>& RHS) const - check for vec-vec and vec-Matrix *|+|-
* vec<Type> Transpose() const - Transpose the vec and return a new vector
* vec<Type>& T() - Transpose the vec and update current vector
* void append(const Type& Ele) - Expand the vec and append the element Ele to the vector
* void resize(const int& Range) - Resize the vec according to the input Range
*********************************************************************************************************************************
*
* To-do:
* Implement an overloaded -= BLAS Level I operator that takes in a RHS vector and update current vector
* Implement a special * operator that inspects the RHS vector, if it is Row * Col, then returns M x N matrix; if it is col * row, then return dot product
*
********************************************************************************************************************************/

#pragma once

#include <boost/multi_array.hpp>
#include <vector>
#include <iostream>

namespace CBLAS
{
    class BLAS
    {
    public:
        BLAS(int M, int N) : m{ M }, n{ N } {} // constructor for vec and 2D matrix
        BLAS(int x, int y, int z) : n1{ x }, n2{ y }, n3{ z } {} // constructor of 3D matrix
        BLAS(int x, int y, int z, int w) : n1{ x }, n2{ y }, n3{ z }, n4{ w } {} // constructor of 4D matrix

    protected:
        int m, n; // M * N size of a 2D matrix/vec
        int n1, n2, n3, n4; // Matrix size for 3D/4D matrix construction
    };

    template <typename Type = double>
    class vec : public BLAS
    {
    public:

        vec() :BLAS{ 0,1 } // Default to 0 row vector
        {
            // V.resize(m);
            dim = "Row";
        }

        explicit vec(int M) : BLAS{ M,1 } // Row based vector constructor
        {
            V.resize(m);
            dim = "Row";
        }

        explicit vec(int M, int N) : BLAS{ M,N }
        {
            if ((m != 1) && (n != 1)) {
                std::cerr << "Invalid instantiation of vec!!" << std::endl;
                exit(-1);
            }
            if ((m == 1) && (n == 0))
            {
                dim = "Col";
            }
            else if ((m == 0) && (n == 1))
            {
                dim = "Row";
            }
            else
            {
                dim = m >= n ? "Row" : "Col";
                int size = m >= n ? m : n;
                V.resize(size);
            }
        }

        constexpr void fill_v(const Type& Val) // fill in the vec with specific value
        {
            for (auto& vec : V) vec = Val;
        }

        Type M() const // return size M
        {
            return m;
        }

        Type N() const // return size N
        {
            return n;
        }

        template <typename A> friend
            std::ostream& operator<<(std::ostream& out, const CBLAS::vec<A>& Vec);

        template<typename A>
        bool check_integrity(const vec<A>& RHS) const // check for vec-vec and vec-Matrix multiplication
        {

            if (*this == RHS)
            {
                type = "Broadcast";
                return true;
            }

            if (m == RHS.N())
            {
                type = "Dot";
                return true;
            }

            return false;

        }

        vec<Type> Transpose() const // Transpose the vec and re-assign
        {
            vec<Type> Vec_T{ n,m };

            Vec_T.V = V;

            return Vec_T;
        }

        vec<Type>& T() // Transpose the intrinsic vec itself
        {
            int m_temp = m;
            m = n;
            n = m_temp;

            if (dim == "Row") dim = "Col";
            else dim = "Row";

            return *this;
        }

        void append(const Type& Ele) // Expand the vec and append the element Ele to the vector 
        {
            if (dim == "Row")
            {
                V.emplace_back(Ele);
                ++m;
            }
            else
            {
                V.emplace_back(Ele);
                ++n;
            }
        }

        void resize(const int& Range) // Resize the vec according to the input Range
        {
            if (dim == "Row")
            {
                m = Range;
                V.resize(m);
            }
            else
            {
                n = Range;
                V.resize(n);
            }
        }

        int size()
        {
            int Dim;
            return Dim = m >= n ? m : n;
        }

        /************************************************************
         *                      Operator Region                     *
        ************************************************************/

        vec<Type>& operator=(const vec<Type>& RHS) // V1 = V - same type
        {
            m = RHS.M();
            n = RHS.N();
            dim = RHS.dim;
            type = RHS.type;

            //V.resize(RHS.V.size());

            V = RHS.V;

            return *this;
        }

        template <typename A>
        vec<Type>& operator=(const vec<A>& RHS) // V1 = V - different type
        {
            m = RHS.M();
            n = RHS.N();
            dim = RHS.dim;
            type = RHS.type;
            V.resize(RHS.V.size());

            for (int i = 0; i < V.size(); ++i)
                V[i] = RHS[i];

            return *this;
        }

        template <typename A>
        vec<Type> operator+(const vec<A>& RHS) const //  V= V1 + V2
        {

            vec<Type> Out{ m, n };

            if (check_integrity<A>(RHS))
            {
                if (type == "Broadcast")
                {
                    for (int i = 0; i < V.size(); ++i)
                        Out[i] = V[i] + RHS[i];

                    return Out;
                }
                else
                {
                    std::cerr << "ERROR:Unmatched input vec dimension for addition!" << std::endl;
                    exit(-1);
                }
            }
            else
            {
                std::cerr << "ERROR:Unmatched input vec dimensions!" << std::endl;
                exit(-1);
            }
        }

        template <typename A>
        vec<Type> operator+(const A const_val) const // V = V1 + a
        {

            vec<Type> Out{ m, n };

            for (int i = 0; i < V.size(); ++i)
                Out[i] = V[i] + const_val;

            return Out;
        }

        template <typename A>
        vec<Type>& operator+=(A const_val) // V1 += a
        {

            for (auto& Vec : V)
                Vec += const_val;

            return *this;
        }

        template <typename A>
        vec<Type>& operator+=(const vec<A>& RHS) // V1 += V regardless whether V1 and V has the same type of dimension
        {

            if (check_integrity<A>(RHS))
            {
                for (int i = 0; i < V.size(); ++i)
                    V[i] += RHS[i];

                return *this;
            }
            else
            {
                std::cerr << "ERROR: Unmatched input size for vector addition!" << std::endl;
                exit(-1);
            }

        }

        template <typename A>
        vec<Type> operator-(const vec<A>& RHS) const //  V= V1 - V2
        {

            vec<Type> Out{ m, n };

            if (check_integrity<A>(RHS))
            {
                if (type == "Broadcast")
                {
                    for (int i = 0; i < V.size(); ++i)
                        Out.V[i] = V[i] - RHS[i];

                    return Out;
                }
                else
                {
                    std::cerr << "ERROR:Unmatched input vec dimension for subtraction!" << std::endl;
                    exit(-1);
                }
            }
            else
            {
                std::cerr << "ERROR:Unmatched input vec dimensions!" << std::endl;
                exit(-1);
            }

        }

        template <typename A>
        vec<Type> operator-(const A const_val) const // V = V1 - a
        {

            vec<Type> Out{ m, n };

            for (int i = 0; i < V.size(); ++i)
                Out[i] = V[i] - const_val;

            return Out;
        }

        template <typename A>
        vec<Type>& operator-=(A const_val) // V1 -= a
        {

            for (auto& Vec : V)
                Vec -= const_val;

            return *this;
        }

        template <typename A>
        vec<Type>& operator-=(const vec<A>& RHS) // V1 -= V regardless whether V1 and V has the same type of dimension
        {
            if (check_integrity<A>(RHS))
            {
                for (int i = 0; i < V.size(); ++i)
                    V[i] -= RHS[i];

                return *this;
            }
            else
            {
                std::cerr << "ERROR: Unmatched input size for vector subtraction!" << std::endl;
                exit(-1);
            }
        }

        template <typename A>
        Type operator*(const vec<A>& RHS) const // Out = V1*V2
        {
            if (check_integrity<A>(RHS))
            {
                Type dot_product{ 0 };
                if (type == "Dot")
                {
                    for (int i = 0; i < V.size(); ++i)
                        dot_product += (V[i] * RHS[i]);
                }
                else
                {
                    std::cout << "Warning: Attempting to perform broadcasting vec multiplication - implicit vec transpose is performed." << std::endl;
                    for (int i = 0; i < V.size(); ++i)
                        dot_product += (V[i] * RHS[i]);
                }
                return dot_product;
            }
            else
            {
                std::cerr << "ERROR:Unmatched input vec dimensions!" << std::endl;
                exit(-1);
            }

        }


        template <typename A>
        vec<Type> operator*(const A const_val) const // V = a * V1
        {
            vec<Type> Out{ m,n };

            for (int i = 0; i < V.size(); ++i)
                Out[i] = V[i] * const_val;

            return Out;
        }

        template <typename A>
        vec<Type>& operator*=(A const_val) // V1 *= a
        {
            for (auto& Vec : V)
                Vec *= const_val;
            return *this;
        }

        Type& operator[](int idx1) // overload the [] operator to return a boost multi-array reference
            // https://theboostcpplibraries.com/boost.multiarray
        {
            return V[idx1];
        }

        const Type& operator[](int idx1) const // overload the [] operator to return a const boost multi-array reference
            // https://theboostcpplibraries.com/boost.multiarray
        {
            return V[idx1];
        }

        /*
        Type& operator()(int idx) // return the ref of value of the idx without knowing the actual layout of the vector
        {
            if (m == 1) return V[0][idx];
            else return V[idx][0];

        }

        const Type& operator()(int idx) const // return the const ref of value of the idx without knowing the actual layout of the vector
        {
            if (m == 1) return V[0][idx];
            else return V[idx][0];
        }
        */

        template <typename A>
        bool operator==(const vec<A>& RHS) const
        {
            if ((m == RHS.M()) && (n == RHS.N())) return true;
            else return false;
        }

        template <typename A>
        bool operator!=(const vec<A>& RHS) const
        {
            if ((m == RHS.M()) && (n == RHS.N())) return false;
            else return true;
        }
    private:
        std::vector<Type> V; // vec is either (m,1) or (1,n)
        std::string dim;            // 1. Row
        // 2. Col
        mutable std::string type;   // 1. Dot - m = n dimension for dot product
        // 2. Broadcast - same m,n dimension for element-wise operation
    };

    template <typename Type>
    std::ostream& operator<<(std::ostream& out, const CBLAS::vec<Type>& Vec)
    {
        if (Vec.dim == "Row")
        {
            for (const auto& Vec_ref : Vec.V)
            {
                out << Vec_ref << std::endl;
            }
        }
        else
        {
            for (const auto& Vec_ref : Vec.V)
                out << Vec_ref << " ";

            out << std::endl;
        }

        return out;
    }

}

template<typename T = double> CBLAS::vec<T> zeros(int m)
{
    CBLAS::vec<T> V{ m,1 };

    return V;
}

template<typename T = double> CBLAS::vec<T> ones(int m)
{
    CBLAS::vec<T> V{ m,1 };
    T val = 1.0;
    V.fill_v(val);
    return V;
}
