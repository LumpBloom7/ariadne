/***************************************************************************
 *            operators.h
 *
 *  Copyright 2008-9  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


/*! \file operators.h
 *  \brief Numerical operator classes
 */

#ifndef ARIADNE_OPERATORS_H
#define ARIADNE_OPERATORS_H

#include <cstdarg>
#include <cassert>
#include <iostream>
#include <string>

#include "tribool.h"
#include "numeric.h"

namespace Ariadne {


enum Comparison {
    EQ,    // Equal
    NEQ,   // Not equal
    GEQ,   // Greater or equal
    LEQ,   // Less or equal
    GT,    // Greater than
    LT,    // Less than
    SGN    // Compare with zero
};

std::ostream& operator<<(std::ostream&, const Comparison& op);

inline std::ostream& operator<<(std::ostream& os, const Comparison& op) {
    switch(op) {
        case SGN:  os << "sgn"; break;
        case EQ:   os << "=="; break;
        case NEQ:  os << "!="; break;
        case GEQ:  os << ">="; break;
        case LEQ:  os << "<="; break;
        case GT:   os << "> "; break;
        case LT:   os << "< "; break;
    }
    return os;
}


enum Operator {
    CNST,  // A constant value
    VAR,   // A named variable
    IND,   // A numbered index
    ADD,   // Addition
    SUB,   // Subtraction
    MUL,   // Multiplication
    DIV,   // Division
    POW,   // Integer power
    POS,   // Unary plus
    NEG,   // Unary negation
    REC,   // Reciprocal
    SQR,   // Square
    SQRT,  // Square root
    EXP,   // Natural exponential
    LOG,   // Natural logarithm
    SIN,   // Sine
    COS,   // Cosine
    TAN,   // Tangent
    ABS,   // Absolute value
    MAX,   // Maximum
    MIN,   // Minimum
    NOT,   // Logical not
    AND,   // Logical and
    OR,    // Logical or
    XOR,   // Logical exclusive or
    IMPL,  // Logical implication
    ITOR   // Conversion of Int to Real
};

std::ostream& operator<<(std::ostream&, const Operator& op);

inline std::ostream& operator<<(std::ostream& os, const Operator& op) {
    switch(op) {
        case CNST: os << "cnst"; break;
        case VAR:  os << "var"; break;
        case IND:  os << "ind"; break;
        case POS:  os << "pos"; break;
        case NEG:  os << "neg"; break;
        case REC:  os << "rec"; break;
        case ADD:  os << "add"; break;
        case SUB:  os << "sub"; break;
        case MUL:  os << "mul"; break;
        case DIV:  os << "div"; break;
        case POW:  os << "pow"; break;
        case NOT:  os << "not"; break;
        case AND:  os << "and"; break;
        case OR:   os << "or"; break;
        case XOR:  os << "xor"; break;
        case IMPL: os << "impl"; break;
        case ABS:  os << "abs"; break;
        case MAX:  os << "max"; break;
        case MIN:  os << "min"; break;
        case SQR:  os << "sqr"; break;
        case SQRT: os << "sqrt"; break;
        case EXP:  os << "exp"; break;
        case LOG:  os << "log"; break;
        case SIN:  os << "sin"; break;
        case COS:  os << "cos"; break;
        case TAN:  os << "tan"; break;
        case ITOR:  os << "itor"; break;
    }
    return os;
}

inline const char* symbol(const Operator& op) {
    switch(op) {
        case POS:  return "+"; break;
        case NEG:  return "-"; break;
        case ADD:  return "+"; break;
        case SUB:  return "-"; break;
        case MUL:  return "*"; break;
        case DIV:  return "/"; break;
        case POW:  return "^"; break;
        case NOT:  return "!"; break;
        case AND:  return "&"; break;
        case OR:   return "|"; break;
        default: assert(false);
    }
}

inline const char* symbol(const Comparison& cmp) {
    switch(cmp) {
        case EQ:  return "=="; break;
        case NEQ: return "!="; break;
        case LEQ: return "<="; break;
        case GEQ: return ">="; break;
        case LT:  return "<"; break;
        case GT:  return ">"; break;
        default: assert(false);
    }
}


struct GtrZero {}; struct LessZero {};

struct Gtr {
    tribool operator()(const Float& x1, const Float& x2) const {
        return (x1==x2) ? indeterminate : tribool(x1>x2); }
    tribool operator()(const Interval& x1, const Interval& x2) const {
        if(x1.lower()>x2.upper()) { return true; } else if(x1.upper()<x2.lower()) { return false; } else { return indeterminate; } }
};

struct Less {
    tribool operator()(const Float& x1, const Float& x2) const {
        return (x1==x2) ? indeterminate : tribool(x1<x2); }
    tribool operator()(const Interval& x1, const Interval& x2) const {
        if(x1.lower()>x2.upper()) { return false; } else if(x1.upper()<x2.lower()) { return true; } else { return indeterminate; } }
};

struct Equal {
    template<class T1, class T2> bool operator()(const T1& a1, const T2& a2) const { return a1 == a2; }
};

struct And { template<class T> T operator()(const T& a1, const T& a2) const { return a1 && a2; } };
struct Or { template<class T> T operator()(const T& a1, const T& a2) const { return a1 || a2; } };
struct Not { template<class T> T operator()(const T& a) const { return !a; } };

struct Add { template<class T> T operator()(const T& a1, const T& a2) const { return a1+a2; } };
struct Sub { template<class T> T operator()(const T& a1, const T& a2) const { return a1-a2; } };
struct Mul { template<class T> T operator()(const T& a1, const T& a2) const { return a1*a2; } };
struct Div { template<class T> T operator()(const T& a1, const T& a2) const { return a1/a2; } };

//struct Pow { template<class T, class N> T operator()(const T& a, const N& n) const { return Ariadne::pow(a,n); } };
struct Pow { Pow(int n) : n(n) { } template<class T> T operator()(const T& a) const { return Ariadne::pow(a,n); } int n; };

struct Neg { template<class T> T operator()(const T& a) const { return Ariadne::neg(a); } };
struct Rec { template<class T> T operator()(const T& a) const { return Ariadne::rec(a); } };
struct Sqr { template<class T> T operator()(const T& a) const { return Ariadne::sqr(a); } };
struct Sqrt { template<class T> T operator()(const T& a) const { return Ariadne::sqrt(a); } };

struct Exp { template<class T> T operator()(const T& a) const { return Ariadne::exp(a); } };
struct Log { template<class T> T operator()(const T& a) const { return Ariadne::log(a); } };
struct Sin { template<class T> T operator()(const T& a) const { return Ariadne::sin(a); } };
struct Cos { template<class T> T operator()(const T& a) const { return Ariadne::cos(a); } };
struct Tan { template<class T> T operator()(const T& a) const { return Ariadne::tan(a); } };


inline std::ostream& operator<<(std::ostream& os, const Less& v) { return os << "<="; }
inline std::ostream& operator<<(std::ostream& os, const Gtr& v) { return os << ">="; }
inline std::ostream& operator<<(std::ostream& os, const Equal& v) { return os << "=="; }

inline std::ostream& operator<<(std::ostream& os, const And& v) { return os << "&&"; }
inline std::ostream& operator<<(std::ostream& os, const Or& v) { return os << "||"; }
inline std::ostream& operator<<(std::ostream& os, const Not& v) { return os << "!"; }

inline std::ostream& operator<<(std::ostream& os, const Add& v) { return os << "+"; }
inline std::ostream& operator<<(std::ostream& os, const Sub& v) { return os << "-"; }
inline std::ostream& operator<<(std::ostream& os, const Mul& v) { return os << "*"; }
inline std::ostream& operator<<(std::ostream& os, const Div& v) { return os << "/"; }

inline std::ostream& operator<<(std::ostream& os, const Neg& op) { return os << "neg"; }
inline std::ostream& operator<<(std::ostream& os, const Rec& op) { return os << "rec"; }
inline std::ostream& operator<<(std::ostream& os, const Pow& op) { return os << "pow"; }
inline std::ostream& operator<<(std::ostream& os, const Sqr& op) { return os << "sqr"; }
inline std::ostream& operator<<(std::ostream& os, const Sqrt& op) { return os << "sqrt"; }

inline std::ostream& operator<<(std::ostream& os, const Exp& op) { return os << "exp"; }
inline std::ostream& operator<<(std::ostream& os, const Log& op) { return os << "log"; }
inline std::ostream& operator<<(std::ostream& os, const Sin& op) { return os << "sin"; }
inline std::ostream& operator<<(std::ostream& os, const Cos& op) { return os << "cos"; }
inline std::ostream& operator<<(std::ostream& os, const Tan& op) { return os << "tan"; }




} // namespace Ariadne

#endif // ARIADNE_OPERATORS_H
