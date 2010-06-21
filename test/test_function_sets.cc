/***************************************************************************
 *            test_function_sets.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include <iostream>

#include "function.h"
#include "box.h"
#include "grid_set.h"
#include "affine_set.h"
#include "function_set.h"
#include "graphics.h"

#include "test.h"

using namespace Ariadne;
using namespace std;


class TestConstrainedImageSet
{
  private:
    Figure figure;
  public:
    void test() {
        figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
        ARIADNE_TEST_CALL(test_draw());
        ARIADNE_TEST_CALL(test_constructor());
        ARIADNE_TEST_CALL(test_geometry());
        ARIADNE_TEST_CALL(test_disjoint());
        ARIADNE_TEST_CALL(test_approximation());
        ARIADNE_TEST_CALL(test_affine_approximation());
        ARIADNE_TEST_CALL(test_split());
        ARIADNE_TEST_CALL(test_draw());
    }

    void test_constructor() {
        List<ScalarFunction> s=ScalarFunction::coordinates(3);
        List<ScalarFunction> x=ScalarFunction::coordinates(2);

        Box d(3,Interval(-1,+2));
        ConstrainedImageSet set(d,(s[0],0.25*s[0]*s[0]+s[1]+0.5*s[2]));
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);
        set.new_space_constraint(x[0]+x[1]<=2.0);
        set.apply((x[0]+x[1],x[0]-x[1]*x[1]));
    }

    void test_geometry() {
        List<ScalarFunction> p=ScalarFunction::coordinates(1);
        List<ScalarFunction> s=ScalarFunction::coordinates(3);
        List<ScalarFunction> x=ScalarFunction::coordinates(2);

        // Test the polytope
        ConstrainedImageSet polytope((Interval(-2,+2),Interval(-2,+2)),(x[0],x[1]));
        polytope.new_parameter_constraint(x[0]+1.5*+x[1]<=1);
        ARIADNE_TEST_ASSERT(polytope.disjoint( (Interval(1.0,2.0),Interval(0.5,1.0)) ));
        ARIADNE_TEST_ASSERT(polytope.overlaps( (Interval(0.0,1.0),Interval(0.5,1.0)) ));

        // Test the unit disc
        ConstrainedImageSet disc((Interval(-2,+2),Interval(-2,+2)),(x[0],x[1]));
        disc.new_parameter_constraint(x[0]*x[0]+x[1]*x[1]<=1);
        ARIADNE_TEST_ASSERT(disc.overlaps( (Interval(-0.5,0.5),Interval(0.25,0.5)) ));
        ARIADNE_TEST_ASSERT(disc.disjoint( (Interval(1,2),Interval(0.5,1)) ));
        ARIADNE_TEST_ASSERT(disc.overlaps( (Interval(0.75,2),Interval(0.5,1)) ));

        // Test a one-dimensional parabolic set
        ConstrainedImageSet parabola(Vector<Interval>(1u,Interval(-1,+1)),(p[0],p[0]*p[0]));
        ARIADNE_TEST_PRINT(parabola);
        ARIADNE_TEST_ASSERT(parabola.disjoint( (Interval(0,0.5),Interval(0.5,1)) ));
        ARIADNE_TEST_ASSERT(parabola.overlaps( (Interval(0.75,2),Interval(0.5,1)) ));

        // Test whether the second iterate of the Henon map intersects a box
        Box d(2,Interval(-0.5,+0.5));
        VectorFunction h((1.5-x[0]*x[0]-0.375*x[1],x[0]));
        VectorFunction f=compose(h,h);
        ConstrainedImageSet set(d,f);
        //set.new_parameter_constraint(0<=x[0]+x[1]<=1);

        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_ASSERT(set.disjoint(Box(2, 1.0,1.25, 1.0,1.25)));
        ARIADNE_TEST_ASSERT(set.overlaps(Box(2, -1.0,-0.875, 1.375,1.625)));
        ARIADNE_TEST_PRINT(f(Point(2,0.375,-0.375)));
    }

    void test_disjoint() {
        List<ScalarFunction> s=ScalarFunction::coordinates(3);
        List<ScalarFunction> x=ScalarFunction::coordinates(2);

        Box d(3,Interval(-1.1,+2.1));
        ConstrainedImageSet set(d,(s[0],0.25*s[0]*s[0]+s[1]+0.5*s[2]));
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);

        Box b(set.bounding_box());
        List<Box> stack(1u,b);
        while(!stack.empty()) {
            Box bx=stack.back();
            stack.pop_back();
            if(bx.radius()>=0.25) {
                Pair<Box,Box> sbx=bx.split();
                stack.append(sbx.first); stack.append(sbx.second);
            } else {
                tribool overlaps = set.overlaps(bx);
                if(overlaps) {
                    figure.set_fill_colour(0.0,0.5,0.5);
                    figure.draw(bx);
                } else if(possibly(overlaps)) {
                    figure.set_fill_colour(0.0,1.0,1.0);
                    figure.draw(bx);
                }
            }
        }
        figure.write("test_constrained_image_set-disjoint");
        figure.clear();
    }

    void test_approximation() {
        List<ScalarFunction> s=ScalarFunction::coordinates(3);
        List<ScalarFunction> x=ScalarFunction::coordinates(2);

        Box d(3,Interval(-1,+2));
        ConstrainedImageSet set(d,(s[0],0.25*s[0]*s[0]+s[1]+0.5*s[2]));
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);
        set.new_space_constraint(x[0]+x[1]<=2.0);
        ARIADNE_TEST_PRINT(set);
        GridTreeSet paving(2);
        uint depth=2;
        paving.adjoin_outer_approximation(set,depth);
        set.adjoin_outer_approximation_to(paving,depth);
        figure.draw(paving);
    }


    void test_split() {
        ScalarFunction o=ScalarFunction::constant(3,1.0);
        ScalarFunction s0=ScalarFunction::coordinate(3,0);
        ScalarFunction s1=ScalarFunction::coordinate(3,1);
        ScalarFunction s2=ScalarFunction::coordinate(3,2);
        ScalarFunction x0=ScalarFunction::coordinate(2,0);
        ScalarFunction x1=ScalarFunction::coordinate(2,1);
        Box d(3,Interval(-1,+1));
        ConstrainedImageSet set(d,(s0,s1+0.5*s2*s2));
        set.new_parameter_constraint(s0+0.75*s1+s2<=0.0);

        ConstrainedImageSet subset1,subset2;
        make_lpair(subset1,subset2)=set.split(0);
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_PRINT(subset1);
        ARIADNE_TEST_PRINT(subset2);
        ARIADNE_TEST_PRINT(set.split(1));

        ConstrainedImageSet subset11,subset12,subset21,subset22;
        make_lpair(subset11,subset12)=subset1.split(0);
        make_lpair(subset21,subset22)=subset2.split(0);
        ARIADNE_TEST_PRINT(subset11);
        subset11.apply(VectorFunction((x0+2.5,x1)));
        ARIADNE_TEST_PRINT(subset11);

        set.apply(VectorFunction((x0-2.5,x1)));
        Figure figure;
        figure.set_bounding_box(Box(2, -4.0,+4.0, -4.0,+4.0));
        figure.set_fill_colour(1.0,1.0,1.0);
        figure.draw(set.bounding_box());
        figure.set_fill_colour(0.75,0.75,0.75);
        figure.draw(set.affine_approximation());
        figure.set_fill_colour(0.5,0.5,0.5);
        figure.draw(subset1.affine_approximation());
        figure.draw(subset2.affine_approximation());
        figure.set_fill_colour(0.25,0.25,0.25);
        figure.draw(subset11.affine_approximation());
        figure.draw(subset12.affine_approximation());
        figure.draw(subset21.affine_approximation());
        figure.draw(subset22.affine_approximation());
        figure.write("test_constrained_image_set-split");
    }

    void test_affine_approximation() {
        // Test conversionn is exact for the affine set -2<x<1; 0<y<2 3x+y<1
        List<ScalarFunction> s=ScalarFunction::coordinates(2);
        Box d(2, -2.0,1.0, 0.0,2.0);
        ConstrainedImageSet set(d,(s[0],s[1]));
        set.new_parameter_constraint(3*s[0]+s[1]<=1);
        AffineSet affine_set=set.affine_approximation();
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_PRINT(affine_set);
        //ARIADNE_TEST_PRINT(set.affine_approximation());
    }

    void test_draw(const std::string& str, const ConstrainedImageSet& set, uint acc) {
        figure.clear();
        GridTreeSet paving(set.dimension());
        set.adjoin_outer_approximation_to(paving,acc);
        figure.set_fill_opacity(1.0);
        figure.set_fill_colour(red);
        figure.draw(paving);
        figure.set_fill_colour(green);
        figure.set_fill_opacity(0.5);
        figure.draw(set);
        figure.write((std::string("test_function_set-draw-")+str).c_str());
        figure.clear();
    }

    void test_draw() {
        ScalarFunction s=ScalarFunction::coordinate(2,0);
        ScalarFunction t=ScalarFunction::coordinate(2,1);
        ScalarFunction x=ScalarFunction::coordinate(2,0);
        ScalarFunction y=ScalarFunction::coordinate(2,1);
        uint acc = 2u;

        test_draw("box",ConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,t)),acc);
        test_draw("polytope",ConstrainedImageSet(Box(2,-2.05,2.05,-1.05,1.05),(s,t),(s+t<=1.5)),acc);
        test_draw("dome",ConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(s,t),(s*s+t<=0.251)),acc);
        test_draw("disc",ConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(s,t),(s*s+t*t<=0.751)),acc);
        test_draw("parallelotope",ConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(2*s+t,s+t)),acc);
        test_draw("ellipse",ConstrainedImageSet(Box(2,-1.0,1.0,-1.0,1.0),(2*s+t,s+t),(s*s+t*t<=0.75)),acc);
        test_draw("concave",ConstrainedImageSet(Box(2,-1.01,1.01,-1.01,1.01),(s,0.25*s*s+t),(2*s+0.25*s*s+t-0.5<=0)),acc);
    }
};




int main(int argc, const char* argv[])
{
    TestConstrainedImageSet().test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}
