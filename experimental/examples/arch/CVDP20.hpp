/***************************************************************************
 *            CVDP20.hpp
 *
 *  Copyright  2020  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "arch.hpp"

using namespace Ariadne;

void CVDP20()
{
    ArchBenchmark benchmark("CVDP20");

    RealVariable x1("x1"), y1("y1"), x2("x2"), y2("y2");

    ARIADNE_LOG_PRINTLN("Coupled van der Pol Oscillator system:")
    ListSet<LabelledEnclosure> reach1, reach2;

    {
        ARIADNE_LOG_PRINTLN_AT(1,"Running for mu=1...");

        RealConstant mu("mu",1.0_dec);
        VectorField dynamics({dot(x1)=y1, dot(y1)=mu*(1-sqr(x1))*y1+x2-2*x1, dot(x2)=y2, dot(y2)=mu*(1-sqr(x2))*y2+x1-2*x2});

        MaximumError max_err = 1e-5;
        TaylorPicardIntegrator integrator(max_err);

        VectorFieldEvolver evolver(dynamics, integrator);
        evolver.configuration().set_maximum_enclosure_radius(0.08);
        evolver.configuration().set_maximum_step_size(0.005);
        evolver.configuration().set_maximum_spacial_error(1e-5);

        RealVariablesBox initial_set({1.25_dec<=x1<=1.55_dec,2.35_dec<=y1<=2.45_dec,1.25_dec<=x2<=1.55_dec,2.35_dec<=y2<=2.45_dec});

        Real evolution_time(7.0);

        StopWatch sw;

        ARIADNE_LOG_PRINTLN_AT(1,"Computing orbit...");
        ARIADNE_LOG_RUN_AT(1,auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER));

        ARIADNE_LOG_PRINTLN_AT(1,"Checking properties...");

        SizeType ce=0;
        for (auto set : orbit.reach()) {
            auto bbox = set.bounding_box();
            if (possibly(bbox[y1] >= 2.75_dec)) {
                ARIADNE_LOG_PRINTLN_AT(2,"set with y1=" << bbox[y1] << " is outside the specification.");
                ++ce;
            }
            if (possibly(bbox[y2] >= 2.75_dec)) {
                ARIADNE_LOG_PRINTLN_AT(2,"set with y2=" << bbox[y2] << " is outside the specification.");
                ++ce;
            }
        }
        sw.click();
        if (ce>0) ARIADNE_LOG_PRINTLN_AT(1,"Number of failures in satisfying the specification: " << ce);
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed() << " seconds.");

        auto instance = benchmark.create_instance("mu1");
        if (ce==0) instance.set_verified(1).set_execution_time(sw.elapsed());
        instance.write();

        reach1.adjoin(orbit.reach());
    }
    {
        ARIADNE_LOG_PRINTLN_AT(1,"Running for mu=2...");

        RealConstant mu("mu",2.0_dec);
        VectorField dynamics({dot(x1)=y1, dot(y1)=mu*(1-sqr(x1))*y1+x2-2*x1, dot(x2)=y2, dot(y2)=mu*(1-sqr(x2))*y2+x1-2*x2});

        MaximumError max_err = 5e-6;
        TaylorPicardIntegrator integrator(max_err);

        VectorFieldEvolver evolver(dynamics, integrator);
        evolver.configuration().set_maximum_enclosure_radius(0.03);
        evolver.configuration().set_maximum_step_size(0.02);
        evolver.configuration().set_maximum_spacial_error(5e-6);

        RealVariablesBox initial_set({1.55_dec<=x1<=1.85_dec,2.35_dec<=y1<=2.45_dec,1.55_dec<=x2<=1.85_dec,2.35_dec<=y2<=2.45_dec});

        Real evolution_time(8.0);

        StopWatch sw;

        ARIADNE_LOG_PRINTLN_AT(1,"Computing orbit...");
        ARIADNE_LOG_RUN_AT(1,auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER));

        ARIADNE_LOG_PRINTLN_AT(1,"Checking properties...");

        SizeType ce=0;
        for (auto set : orbit.reach()) {
            auto bbox = set.bounding_box();
            if (possibly(bbox[y1] >= 4.05_dec)) {
                ARIADNE_LOG_PRINTLN_AT(2,"set with y1=" << bbox[y1] << " is outside the specification.");
                ++ce;
            }
            if (possibly(bbox[y2] >= 4.05_dec)) {
                ARIADNE_LOG_PRINTLN_AT(2,"set with y2=" << bbox[y2] << " is outside the specification.");
                ++ce;
            }
        }
        sw.click();
        if (ce>0) ARIADNE_LOG_PRINTLN_AT(1,"Number of failures in satisfying the specification: " << ce);
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed() << " seconds.");

        auto instance = benchmark.create_instance("mu2");
        if (ce==0) instance.set_verified(1).set_execution_time(sw.elapsed());
        instance.write();

        reach2.adjoin(orbit.reach());
    }

    ARIADNE_LOG_PRINTLN("Plotting...");

    LabelledFigure fig(Axes2d(-2.5<=x1<=2.5,-4.05<=y1<=4.05));
    fig << fill_colour(ariadneorange);
    fig.draw(reach2);
    fig << fill_colour(Colour(0.6,0.6,0.6));
    fig.draw(reach1);
    fig.write(benchmark.name().c_str());
    ARIADNE_LOG_PRINTLN("File " << benchmark.name() << ".png written.");
}
