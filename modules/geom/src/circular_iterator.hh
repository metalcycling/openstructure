//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------
#ifndef CIRCULAR_ITERATOR_HH
#define CIRCULAR_ITERATOR_HH

#include <iterator>

template<typename T>
class const_circular_iter : public std::iterator<std::bidirectional_iterator_tag,typename T::value_type>  
{
typedef typename T::const_iterator iterator;
typedef typename T::value_type value_type;
protected:
        iterator iter;
        iterator begin;
        iterator end;
public:
    const_circular_iter(T& t) : iter(t.begin()), begin(t.begin()), end(t.end()) {};
    const_circular_iter(iterator b,iterator e) : iter(b), begin(b), end(e) {};
    const_circular_iter(iterator b,iterator e,iterator pos) : iter(pos), begin(b), end(e) {};
    operator iterator(){
      return iter;
    }
    const_circular_iter<T>& operator++()
    {   
        ++iter;

        if (iter == end)
            iter = begin;

        return(*this);
    }

    const_circular_iter<T> operator++(int)
    {
        const_circular_iter<T>  t=*this;
        ++(*this);
        return(t);
    }
    const_circular_iter<T>& operator--()
    {   
        if (iter == begin)
            iter = end;

        --iter;

        return(*this);
    }

    const_circular_iter<T> operator--(int)
    {
        const_circular_iter<T>  t=*this;
        --(*this);
        return(t);
    }

    const_circular_iter<T> operator-(int n)
    {
        const_circular_iter<T>  t=*this;
        for(unsigned int i=0;i<n;++i)
        --t;
        return t;
    }
    const_circular_iter<T> operator+(int n)
    {
        const_circular_iter<T>  t=*this;
        for(unsigned int i=0;i<n;++i)
        ++t;
        return t;
    }

    const value_type operator*() const 
        {return (*iter);}

    bool operator==(const const_circular_iter<T>&   rhs) const
        {return (iter == rhs.iter);}

    bool operator!=(const const_circular_iter<T>&   rhs) const 
        {return ! operator==(rhs); }
};


template<typename T>
class circular_iter : public std::iterator<std::bidirectional_iterator_tag,typename T::value_type>  
{
typedef typename T::iterator iterator;
typedef typename T::value_type value_type;
private:
        iterator     iter;
        iterator     begin;
        iterator     end;
public:
    circular_iter(T& t) : iter(t.begin()), begin(t.begin()), end(t.end()) {};
    circular_iter(iterator b, iterator e) : iter(b), begin(b), end(e) {};
    circular_iter(iterator b, iterator e, iterator pos) : iter(pos), begin(b), end(e) {};

    circular_iter<T>& operator++()
    {   
        ++iter;

        if (iter == end)
            iter = begin;

        return(*this);
    }

    circular_iter<T> operator++(int)
    {
        circular_iter<T>    t=*this;
        ++(*this);
        return(t);
    }
    circular_iter<T>& operator--()
    {   
        if (iter == begin)
            iter = end;

        --iter;

        return(*this);
    }

    circular_iter<T> operator--(int)
    {
        circular_iter<T>    t=*this;
        --(*this);
        return(t);
    }
    circular_iter<T> operator-(int n)
    {
      circular_iter<T>  t=*this;
        for(unsigned int i=0;i<n;++i)
        --t;
        return t;
    }
    circular_iter<T> operator+(int n)
    {
      circular_iter<T>  t=*this;
        for(unsigned int i=0;i<n;++i)
        ++t;
        return t;
    }

    const value_type operator*() const 
    {return (*iter);}

    value_type& operator*() 
    {return (*iter);}

    bool operator==(const circular_iter<T>& rhs) const
    {return (iter == rhs.iter);}

    bool operator!=(const circular_iter<T>& rhs) const {return ! operator==(rhs); }
    operator iterator(){
      return iter;
    }
};

#endif

