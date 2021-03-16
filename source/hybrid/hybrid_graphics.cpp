/***************************************************************************
 *            hybrid/hybrid_graphics.cpp
 *
 *  Copyright  2011-20  Pieter Collins
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

#include "function/functional.hpp"
#include "config.hpp"

#include "utility/macros.hpp"
#include "utility/stlio.hpp"
#include "numeric/numeric.hpp"
#include "symbolic/space.hpp"
#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "output/geometry2d.hpp"
#include "output/graphics.hpp"
#include "hybrid/discrete_location.hpp"
#include "geometry/function_set.hpp"
#include "symbolic/expression_set.hpp"
#include "hybrid/hybrid_graphics.hpp"

namespace Ariadne {

static const Nat HYBRID_DEFAULT_WIDTH = 800;
static const Nat HYBRID_DEFAULT_HEIGHT = 800;

static const Nat HYBRID_LEFT_MARGIN = 0;
static const Nat HYBRID_BOTTOM_MARGIN = 0;
static const Nat HYBRID_TOP_MARGIN = 0;
static const Nat HYBRID_RIGHT_MARGIN = 0;


HybridFigure::~HybridFigure() {
}


HybridFigure::HybridFigure()
    : variables(RealVariable("x"),RealVariable("y"))
{
}

Void HybridFigure::write(const char* cfilename) const
{
    #ifdef HAVE_CAIRO_H
        this->write(cfilename, CairoFileType::PNG);
    #else
    #ifdef HAVE_GNUPLOT_H
        this->write(cfilename, GnuplotFileType::PNG);
    #else
        ARIADNE_ERROR("No facilities for displaying graphics are available.");
    #endif
    #endif
}

Void
HybridFigure::write(const char* cfilename, CairoFileType fileType) const
{
    this->write(cfilename, HYBRID_DEFAULT_WIDTH, HYBRID_DEFAULT_HEIGHT, fileType);
}

Void
HybridFigure::write(const char* cfilename, GnuplotFileType fileType) const
{
    this->write(cfilename, HYBRID_DEFAULT_WIDTH, HYBRID_DEFAULT_HEIGHT, fileType);
}


Void
HybridFigure::write(const char* cfilename, Nat drawing_width, Nat drawing_height, CairoFileType fileType) const
{
    const Nat canvas_width = drawing_width+HYBRID_LEFT_MARGIN+HYBRID_RIGHT_MARGIN;
    const Nat canvas_height = drawing_height+HYBRID_BOTTOM_MARGIN+HYBRID_TOP_MARGIN;

    SharedPointer<CanvasInterface> canvas=make_canvas(cfilename, canvas_width, canvas_height, fileType);

    this->_paint_all(*canvas);
/*
    StringType filename(cfilename);
    if(filename.rfind(".") != StringType::npos) {
    } else {
        filename=filename+".png";
    }
*/
    canvas->write(cfilename);
}

Void
HybridFigure::write(const char* cfilename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType) const
{
    const Nat canvas_width = drawing_width+HYBRID_LEFT_MARGIN+HYBRID_RIGHT_MARGIN;
    const Nat canvas_height = drawing_height+HYBRID_BOTTOM_MARGIN+HYBRID_TOP_MARGIN;

    SharedPointer<CanvasInterface> canvas=make_canvas(cfilename, canvas_width, canvas_height, fileType);

    this->_paint_all(*canvas);
/*
    StringType filename(cfilename);
    if(filename.rfind(".") != StringType::npos) {
    } else {
        filename=filename+".png";
    }
*/
    canvas->write(cfilename);
}


Void HybridFigure::_paint_all(CanvasInterface& canvas) const
{
    // Project the bounding box onto the canvas
    double xl=numeric_cast<double>(bounds[variables.x_variable()].lower_bound());
    double xu=numeric_cast<double>(bounds[variables.x_variable()].upper_bound());
    double yl=numeric_cast<double>(bounds[variables.y_variable()].lower_bound());
    double yu=numeric_cast<double>(bounds[variables.y_variable()].upper_bound());

    canvas.initialise(variables.x_variable().name(),variables.y_variable().name(),xl,xu,yl,yu);

    // Draw shapes
    for(SizeType i=0; i!=objects.size(); ++i) {
        const HybridDrawableInterface& shape=*objects[i].shape_ptr;
        set_properties(canvas, objects[i].properties);
        shape.draw(canvas,this->locations,this->variables);
    }

    canvas.finalise();
}



} // namespace Ariadne


