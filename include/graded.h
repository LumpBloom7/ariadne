/***************************************************************************
 *            graded.h
 *
 *  Copyright  2011  Pieter Collins
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

/*! \file graded.h
 *  \brief Graded algebras.
 */

#ifndef ARIADNE_GRADED_H
#define ARIADNE_GRADED_H

#include "procedure.h"

namespace Ariadne {

struct AntiDiff { };

template<class Op, class A1, class A2=Void, class A3=Void> struct ClosureExpression;

template<class Op, class A> struct ClosureExpression<Op,A,Void,Void> {
    const A& arg;
    ClosureExpression(const A& a) : arg(a) { }
};

template<class Op, class A1, class A2> struct ClosureExpression<Op,A1,A2,Void> {
    const A1& arg1; const A2& arg2;
    ClosureExpression(const A1& a1, const A2& a2) : arg1(a1), arg2(a2) { }
};

template<class Op, class A1, class A2, class A3> struct ClosureExpression {
    const A1& arg1; const A2& arg2; const A2& arg3;
    ClosureExpression(const A1& a1, const A2& a2, const A3& a3) : arg1(a1), arg2(a2), arg3(a3) { }
};

template<class Op, class A> ClosureExpression<Op,A>
make_expression(Op op, const A& a) {
    return ClosureExpression<Op,A>(a);
};

template<class Op, class A1, class A2> ClosureExpression<Op,A1,A2>
make_expression(Op op, const A1& a1, const A2& a2) {
    return ClosureExpression<Op,A1,A2>(a1,a2);
};

template<class Op, class A1, class A2, class A3> ClosureExpression<Op,A1,A2,A3>
make_expression(Op op, const A1& a1, const A2& a2, const A3& a3) {
    return ClosureExpression<Op,A1,A2,A3>(a1,a2,a3);
};



bool compatible(const Float& x1, const Float& x2) { return true; }
template<class X> bool compatible(const Polynomial<X>& x1, const Polynomial<X>& x2) { return x1.argument_size()==x2.argument_size(); }
template<class X> bool compatible(const Differential<X>& x1, const Differential<X>& x2) { return x1.argument_size()==x2.argument_size(); }

Float create(const Float& x) { return Float(0); }
Interval create(const Interval& x) { return Interval(0); }
template<class X> Polynomial<X> create(const Polynomial<X>& x) { return Polynomial<X>(x.argument_size()); }
template<class X> Differential<X> create(const Differential<X>& x) { return Differential<X>(x.argument_size(),x.degree()); }

std::ostream& operator<<(std::ostream& os, const Differential<Interval>& d) {
    return os << Polynomial<Interval>(d.expansion());
}

template<class A> struct Graded : public List<A>
{
    typedef Graded<A> SelfType;
    Graded() : List<A>() { }
    Graded(const A& a) : List<A>(1u,a) { }
    Graded(const List<A>& l) : List<A>(l) { }
    template<class Op> Void operator=(const ClosureExpression<Op,SelfType,SelfType>& expr);
    template<class Op> Void operator=(const ClosureExpression<Op,SelfType>& expr);
    Void operator=(const ClosureExpression<AntiDiff,SelfType>& ad);
    Nat degree() const { return this->size()-1u; }
    Void extend(const A& a) { this->List<A>::append(a); }
};
template<class A> std::ostream& operator<<(std::ostream& os, const Graded<A>& g) {
    if(g.size()==0) { return os << "0"; }
    os << "(" << g[0] << ")";
    for(uint i=1; i<=g.degree(); ++i) {
        os << " + (" << g[i] << ")*t";
        if(i>1) { os << "^"<<i; }
    }
    return os;
}
template<> std::ostream& operator<<(std::ostream& os, const Graded<Float>& g) {
    if(g.size()==0) { return os << "0"; }
    os << g[0];
    for(uint i=1; i<=g.degree(); ++i) {
        if(g[i]>0) {
            if(g[i]==1) { os << "+t"; } else { os << "+"<<g[i]<<"*t"; }
            if(i>1) { os << "^"<<i; }
        } else if(g[i]<0) {
            if(g[i]==-1) { os << "-t"; } else { os << g[i]<<"*t"; }
            if(i>1) { os << "^"<<i; }
        }
    }
    return os;
}
template<> std::ostream& operator<<(std::ostream& os, const Graded<Interval>& g) {
    if(g.size()==0) { return os << "0"; }
    os << g[0];
    for(uint i=1; i<=g.degree(); ++i) {
        os << "+"<<g[i]<<"*t";
        if(i>1) { os << "^"<<i; }
    }
    return os;
}

template<class X> Void compute(X& r, const Add&, const X& a1, const X& a2) { return add(r,a1,a2); }
template<class X> Void compute(X& r, const Sub&, const X& a1, const X& a2) { return sub(r,a1,a2); }
template<class X> Void compute(X& r, const Mul&, const X& a1, const X& a2) { return mul(r,a1,a2); }
template<class X> Void compute(X& r, const Div&, const X& a1, const X& a2) { return div(r,a1,a2); }
template<class X> Void compute(X& r, const Sqr&, const X& a) { return sqr(r,a); }
template<class X> Void compute(X& r, const Sqrt&, const X& a) { return sqrt(r,a); }
template<class X> Void compute(X& r, const Exp&, const X& a) { return exp(r,a); }
template<class X> Void compute(X& r, const Log&, const X& a) { return log(r,a); }
template<class X> Void compute(X& r, const Rec&, const X& a) { return rec(r,a); }


template<class A, class B> Graded<A>& operator+=(Graded<A>& a, const B& c) {
    if(a.degree()==1) { a[0]+=c; } return a; }
template<class A, class B> Graded<A>& operator*=(Graded<A>& a, const B& c) {
    a.back()*=c; return a; }

template<class A, class B> Graded<A> operator+(const Graded<A>& a, const B& c) {
    Graded<A> r(a); r+=c; return r;  }
template<class A> Graded<A> operator*(const Graded<A>& a, int c) {
    Graded<A> r(a); r*=c; return r;  }

template<class A> ClosureExpression<Neg,Graded<A> > operator-(const Graded<A>& a) {
    return make_expression(Neg(),a); }
template<class A> ClosureExpression<Add,Graded<A>,Graded<A> > operator+(const Graded<A>& a1, const Graded<A>& a2) {
    return make_expression(Add(),a1,a2); }
template<class A> ClosureExpression<Sub,Graded<A>,Graded<A> > operator-(const Graded<A>& a1, const Graded<A>& a2) {
    return make_expression(Sub(),a1,a2); }
template<class A> ClosureExpression<Mul,Graded<A>,Graded<A> > operator*(const Graded<A>& a1, const Graded<A>& a2) {
    return make_expression(Mul(),a1,a2); }
template<class A> ClosureExpression<Div,Graded<A>,Graded<A> > operator/(const Graded<A>& a1, const Graded<A>& a2) {
    return make_expression(Div(),a1,a2); }
template<class A> ClosureExpression<Sqr,Graded<A> > sqr(const Graded<A>& a) {
    return make_expression(Sqr(),a); }
template<class A> ClosureExpression<Sqrt,Graded<A> > sqrt(const Graded<A>& a) {
    return make_expression(Sqrt(),a); }
template<class A> ClosureExpression<Exp,Graded<A> > exp(const Graded<A>& a) {
    return make_expression(Exp(),a); }
template<class A> ClosureExpression<Log,Graded<A> > log(const Graded<A>& a) {
    return make_expression(Log(),a); }
template<class A> ClosureExpression<Rec,Graded<A> > rec(const Graded<A>& a) {
    return make_expression(Rec(),a); }

template<class A> Graded<A> pos(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> neg(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> abs(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> pow(const Graded<A>& a, int n) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> sin(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> cos(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> tan(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> atan(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }

template<class A> Void add(Graded<A>& r, const Graded<A>& a1, const Graded<A>& a2) {
    ARIADNE_ASSERT(r.degree()+1u == a1.degree());
    ARIADNE_ASSERT(r.degree()+1u == a2.degree());
    r.append(a1.back()+a2.back());
}

template<class A> Void sub(Graded<A>& r, const Graded<A>& a1, const Graded<A>& a2) {
    ARIADNE_ASSERT(r.degree()+1u == a1.degree());
    ARIADNE_ASSERT(r.degree()+1u == a2.degree());
    r.append(a1.back()-a2.back());
}

template<class A> Void mul(Graded<A>& r, const Graded<A>& a1, const Graded<A>& a2) {
    ARIADNE_ASSERT(r.degree()+1u == a1.degree());
    ARIADNE_ASSERT(r.degree()+1u == a2.degree());
    ARIADNE_ASSERT(compatible(a1[0],a2[0]));
    r.append(create(a1[0]));
    Nat d = r.degree();
    for(Nat i=0; i<=d; ++i) {
        r[d] += a1[i]*a2[d-i];
    }
}

template<class A> Void div(Graded<A>& r, const Graded<A>& a1, const Graded<A>& a2) {
    ARIADNE_ASSERT(r.degree()+1u == a1.degree());
    ARIADNE_ASSERT(r.degree()+1u == a2.degree());
    ARIADNE_ASSERT(compatible(a1[0],a2[0]));
    r.append(create(a1[0]));
    Nat d = r.degree();
    r[d]+=a1[d];
    for(Nat i=0; i!=d; ++i) {
        r[d] -= a2[d-i]*r[i];
    }
    r[d]=r[d]/a2[0];
}

template<class A> Void sqr(Graded<A>& r, const Graded<A>& a) {
    std::cerr<<"r="<<r<<" a="<<a<<"\n";
    ARIADNE_ASSERT(r.size() < a.size());
    r.append(create(a[0]));
    Nat d = r.degree();
    for(Nat i=0; i<=d; ++i) {
        r[d] += a[i]*a[d-i];
    }
}

template<class A> Void rec(Graded<A>& r, const Graded<A>& a) {
    ARIADNE_ASSERT(r.degree()+1u == a.degree());
    r.append(create(a[0]));
    Nat d = r.degree();
    if(d==0) { r[d]=rec(a[0]); return; }
    for(Nat i=0; i!=d; ++i) {
        r[d] -= a[d-i]*r[i];
    }
    r[d]=r[d]*r[0];
}

template<class A> Void sqrt(Graded<A>& r, const Graded<A>& a) {
    ARIADNE_ASSERT(r.degree()+1u == a.degree());
    r.append(create(a[0]));
    Nat d = r.degree();
    if(d==0) { r[d]=sqrt(a[0]); return; }
    r[d]=a[d];
    for(Nat i=1; i!=d; ++i) {
        r[d] -= r[d-i]*r[i];
    }
    r[d]=r[d]/(2*r[0]);
}

template<class A> Void exp(Graded<A>& r, const Graded<A>& a) {
    ARIADNE_ASSERT(r.degree()+1u == a.degree());
    r.append(create(a[0]));
    Nat d = r.degree();
    if(d==0) { r[d]+=exp(a[0]); return; }
    for(Nat i=0; i!=d; ++i) {
        r[d] += (d-i)*a[d-i]*r[i];
    }
    r[d]/=d;
}

template<class A> Void log(Graded<A>& r, const Graded<A>& a) {
    // y=log x; r=1/x; s=r^2;
    r.append(create(a[0]));
    Nat d = r.degree();
    if(d==0) { r[d]=log(a[0]); return; }
    if(d==1) { r[d]=a[d]/a[0];  return; }
    r[d]+=d*a[d];
    for(Nat i=1; i!=d; ++i) {
        r[d] -= a[d-i]*r[i]*i;
    }
    r[d]=r[d]/a[0];
    r[d]/=d;
}


template<class A> template<class Op> Void Graded<A>::operator=(const ClosureExpression<Op,Graded<A>,Graded<A> >& expr) {
    compute(*this,Op(),expr.arg1,expr.arg2);
}

template<class A> template<class Op> Void Graded<A>::operator=(const ClosureExpression<Op,Graded<A> >& expr) {
    compute(*this,Op(),expr.arg);
}

template<class A> Void Graded<A>::operator=(const ClosureExpression<AntiDiff,Graded<A> >& expr) {
    antidifferentiate(*this,expr.arg);
}

template<class A> Void antidifferentiate(Graded<A>& r, const Graded<A>& a) {
    ARIADNE_ASSERT(r.degree() == a.degree());
    r.append(a.back()/(r.degree()+1u));
}

template<class A> ClosureExpression<AntiDiff,Graded<A> > antidifferential(const Graded<A>& a) {
    return ClosureExpression<AntiDiff,Graded<A> >(a);
}


template<class T>
Vector< Polynomial<T> > flow(const Vector< Procedure<T> >& vf, Vector<T>& dom)
{
    List< Polynomial<T> > intermediates(vf._instructions.size(), Polynomial<T>(dom.size()));
}

Pair<List<Float>,Float> midpoint_error(const Graded<Interval>& x) {
    List<Float> m(x.degree()+1);
    Float e;
    for(uint i=0; i<=x.degree(); ++i) {
        m[i]=midpoint(x[i]);
        e=add_up(e,max(sub_up(m[i],x[i].lower()),sub_up(x[i].upper(),m[i])));
    }
    return Pair<List<Float>,Float>(m,e);
}


template<class X> Differential<X> make_differential_variable(Nat as, Nat deg, X val, Nat ind) {
    return Differential<X>::variable(as,deg,val,ind); }
template<class X> Graded<X> make_graded(const X& val) {
    return Graded<X>(val); }
template<class X> Graded<X> create_graded(const X&) {
    return Graded<X>(); }

struct Function {
    template<class X> X operator()(const X& x) const { return x; }
};

Float eval(const Graded<Float>& p, const Float& x) {
    Nat d=p.degree();
    Float r=p[d];
    for(uint i=0; i!=d; ++i) {
        r=r*x+p[d-i-1];
    }
    std::cerr<<"eval("<<p<<","<<x<<")="<<r<<"\n";
    return r;
}



template<class A> void compute(const Vector< Procedure<Real> >& p, Vector< Graded<A> >& r, List< Graded<A> >& t, const Vector< Graded<A> >& a) {
    _compute(t,p._instructions,p._constants,a);
    for(uint i=0; i!=p._results.size(); ++i) { r[i]=t[p._results[i]]; }
}

template<class A> Vector< Graded<A> > integrate(const Vector< Procedure<Real> >& p, const Vector<A>& x) {
    Vector< Graded<A> > arg(x.size(),create(x[0]));
    for(uint i=0; i!=x.size(); ++i) { arg[i]=Graded<A>(x[i]); }
    std::cerr<< "x0="<<arg <<"\n";

    List< Graded<A> > tmp(p.temporaries_size());
     Vector< Graded<A> > res(p.result_size(),create(x[0]));

    const uint N = 6;

    for(uint n=0; n!=N; ++n) {
        std::cerr<<"arg="<<arg<<"\n";
        compute(p,res,tmp,arg);
        std::cerr<<"    tmp="<<tmp<<"\n";
        std::cerr<<"  res="<<res<<"\n";
        for(uint i=0; i!=arg.size(); ++i) {
            arg[i]=antidifferential(res[i]);
        }
    }
    std::cerr<<"arg="<<arg<<"\n";

    return arg;
}

Vector< Graded< Differential<Interval> > > integrate(const Vector< Procedure<Real> >& p, const Vector<Interval>& x) {
    const uint M = 2;
    const uint N = 6;
    Vector< Graded<Differential<Interval> > > arg(x.size(),Differential<Interval>(x.size(),M));
    //for(uint i=0; i!=x.size(); ++i) { arg[i]=Differential<Interval>::variable(x.size(),M,x[i],i); }
    for(uint i=0; i!=x.size(); ++i) { arg[i]=Differential<Interval>::constant(x.size(),M,x[i]); }

    Graded< Differential<Interval> > null;

    List< Graded<Differential<Interval> > > tmp(p.temporaries_size(),null);
    Vector< Graded<Differential<Interval> > > res(p.result_size(),null);

    std::cerr<< "arg="<<arg <<"\n";
    std::cerr<< "tmp="<<tmp <<"\n";
    std::cerr<< "res="<<res <<"\n";

    for(uint n=0; n!=N; ++n) {
        std::cerr<<"arg="<<arg<<"\n";
        compute(p,res,tmp,arg);
        std::cerr<<"    tmp="<<tmp<<"\n";
        std::cerr<<"  res="<<res<<"\n";
        for(uint i=0; i!=arg.size(); ++i) {
            arg[i]=antidifferential(res[i]);
        }
    }
    std::cerr<<"arg="<<arg<<"\n";

    return arg;
}


Vector< Graded< Differential<Interval> > > flow(const Vector< Procedure<Real> >& p, const Vector<Interval>& x) {
    const uint M = 2;
    const uint N = 6;

    const uint AS = x.size();

    IntervalVector D=x;
    IntervalVector C=midpoint(x);

    //IntervalVector B=D;
    //IntervalVector A=C;

    Graded< Differential<Interval> > null;
    List< Graded<Differential<Interval> > > tmp(p.temporaries_size(),null);
    Vector< Graded<Differential<Interval> > > res(p.result_size(),null);

    Vector< Graded<Differential<Interval> > > arg(x.size(),Differential<Interval>(x.size(),M));
    for(uint i=0; i!=x.size(); ++i) { arg[i]=Differential<Interval>::variable(x.size(),M,D[i],i); }


    std::cerr<< "arg="<<arg <<"\n";
    std::cerr<< "tmp="<<tmp <<"\n";
    std::cerr<< "res="<<res <<"\n";

    for(uint n=0; n!=N; ++n) {
        std::cerr<<"arg="<<arg<<"\n";
        compute(p,res,tmp,arg);
        std::cerr<<"    tmp="<<tmp<<"\n";
        std::cerr<<"  res="<<res<<"\n";
        for(uint i=0; i!=arg.size(); ++i) {
            arg[i]=antidifferential(res[i]);
        }
    }
    std::cerr<<"arg="<<arg<<"\n";

    return arg;
}

} // namespace Ariadne

#endif /* ARIADNE_GRADED_H */