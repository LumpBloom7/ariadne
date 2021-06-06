/***************************************************************************
 *            noisy-utilities.hpp
 *
 *  Copyright  2008-18 Luca Geretti
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

#include "ariadne.hpp"
#include "utility/stopwatch.hpp"

namespace Ariadne {

typedef Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> SystemType;
typedef Real ContinuousTimeType;

void run_single(String name, DifferentialInclusion const& ivf, BoxDomainType const& initial, ContinuousTimeType evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, Reconditioner const& reconditioner);
void run_each_approximation(String name, DifferentialInclusion const& ivf, BoxDomainType const& initial, ContinuousTimeType evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, Reconditioner const& reconditioner);

void run_noisy_system(String name, DottedRealAssignments const& dynamics, RealVariablesBox const& inputs, RealVariablesBox const& initial, ContinuousTimeType evolution_time, double step);
void run_noisy_system(SystemType system);


template<class F, class S> List<ResultOf<F(S)>> map(F const& f, List<S> const& list) {
    List<ResultOf<F(S)>> result; for(auto item : list) { result.append(f(item)); } return result;
}

inline ApproximateDouble score(ValidatedConstrainedImageSet const& evolve_set) {
    auto bbx = evolve_set.bounding_box();
    return 1.0/std::pow(volume(bbx).get_d(),1.0/bbx.size());
}

void run_single(String name, DifferentialInclusion const& ivf, BoxDomainType const& initial, Real evolution_time, ApproximateDouble step, List<InputApproximation> approximations, SweeperDP sweeper, IntegratorInterface const& integrator, Reconditioner const& reconditioner) {
    CONCLOG_SCOPE_CREATE;
    auto evolver = DifferentialInclusionEvolver(ivf, sweeper, integrator, reconditioner);
    evolver.configuration().approximations(approximations);
    evolver.configuration().maximum_step_size(step);

    Stopwatch<Milliseconds> sw;

    List<ValidatedVectorMultivariateFunctionPatch> flow_functions = evolver.reach(initial,evolution_time);
    sw.click();

    List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorMultivariateFunctionPatch const& fm){return ValidatedConstrainedImageSet(fm.domain(),fm);},flow_functions);
    auto final_set = flow_functions.back();
    ValidatedVectorMultivariateFunctionPatch evolve_function =
        partial_evaluate(final_set,final_set.result_size(),static_cast<ExactNumber>(final_set.domain()[final_set.result_size()].upper_bound()));
    auto evolve_set = ValidatedConstrainedImageSet(evolve_function.domain(),evolve_function);

    CONCLOG_PRINTLN("Score: " << score(evolve_set) << ", time: " << sw.elapsed_seconds() << " s");

    CONCLOG_PRINTLN("Plotting...");
    auto n = ivf.dimension();
    Box<FloatDPUpperInterval> graphics_box(n);
    for (auto set: reach_sets) {
        graphics_box = hull(graphics_box,set.bounding_box());
    }
    for (SizeType i : range(0,n-1)) {
        for (SizeType j : range(i+1,n)) {
            Figure fig=Figure(graphics_box,Projection2d(n,i,j));
            for (auto set : reach_sets) { fig.draw(set); }
            char num_char[64] = "";
            if (n > 2) sprintf(num_char,"[%lu,%lu]",i,j);
            CONCLOG_RUN_AT(1,fig.write((name+num_char).c_str()));
        }
    }
}

void run_each_approximation(String name, DifferentialInclusion const& ivf, BoxDomainType const& initial, Real evolution_time, double step, List<InputApproximation> approximations, SweeperDP sweeper, IntegratorInterface const& integrator, Reconditioner const& reconditioner) {
    for (auto appro: approximations) {
        List<InputApproximation> singleapproximation = {appro};
        CONCLOG_PRINTLN(appro);
        run_single(name,ivf,initial,evolution_time,step,singleapproximation,sweeper,integrator,reconditioner);
    }
}

void run_noisy_system(String name, const DottedRealAssignments& dynamics, const RealVariablesBox& inputs,
              const RealVariablesBox& initial, Real evolution_time, double step) {

    DifferentialInclusion ivf(dynamics, inputs);

    SizeType period_of_parameter_reduction=12;
    ExactDouble ratio_of_parameters_to_keep=6.0_x;
    double sw_threshold = 1e-8;
    ThresholdSweeperDP sweeper(DoublePrecision(),sw_threshold);

    List<InputApproximation> approximations;
    approximations.append(ZeroApproximation());
    approximations.append(ConstantApproximation());
    approximations.append(AffineApproximation());
    approximations.append(SinusoidalApproximation());
    approximations.append(PiecewiseApproximation());

    TaylorPicardIntegrator integrator(
            step_maximum_error=1e-3,
            sweeper,
            lipschitz_tolerance=0.5_x,
            minimum_temporal_order=4,
            maximum_temporal_order=12);

    LohnerReconditioner reconditioner(initial.variables().size(),inputs.variables().size(),period_of_parameter_reduction,ratio_of_parameters_to_keep);

    run_single(name,ivf,initial_ranges_to_box(initial),evolution_time,step,approximations,sweeper,integrator,reconditioner);
    //run_each_approximation(name,ivf,initial_ranges_to_box(initial),evolution_time,step,approximations,sweeper,integrator,reconditioner);
}

void run_noisy_system(SystemType system) {
    run_noisy_system(std::get<0>(system),std::get<1>(system),std::get<2>(system),std::get<3>(system),std::get<4>(system),std::get<5>(system));
}

}
