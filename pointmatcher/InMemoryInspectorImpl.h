// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Copyright (c) 2010--2012,
François Pomerleau and Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the authors at <f dot pomerleau at gmail dot com> and
<stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef __POINTMATCHER_IN_MEMORY_INSPECTORS_H
#define __POINTMATCHER_IN_MEMORY_INSPECTORS_H

#include "PointMatcher.h"

namespace PointMatcherSupport
{

/**
 * Why is InMemoryInspector not in registry and not in InspectorImpl?
 * Contrary to other module of libpointmatcher, InMemoryInspector is a module that should only
 * be config into Libpm programmatically.  The goal of this inspector is to log all the stats
 * and copy them into memory. Then the user can read those stats direcly into memory, this
 * remove the requirement of dump the stats to the disk and then parsing the dump.
 *
 * The user must have access to the implementation of InMemoryInspector otherwise it can not
 * read it's member function and get the stats. In libpointmatcher the implementation of other
 * module's are not public to the user.
 */
template<typename T>
struct InMemoryInspector: public PointMatcher<T>::Inspector
{
	 typedef PointMatcherSupport::Parametrizable Parametrizable;
	 typedef PointMatcherSupport::Parametrizable P;
	 typedef Parametrizable::Parameters Parameters;
	 typedef Parametrizable::ParameterDoc ParameterDoc;
	 typedef Parametrizable::ParametersDoc ParametersDoc;

	 typedef typename PointMatcher<T>::Inspector Inspector;
	 typedef typename PointMatcher<T>::DataPoints DataPoints;
	 typedef typename PointMatcher<T>::Matches Matches;
	 typedef typename PointMatcher<T>::OutlierWeights OutlierWeights;
	 typedef typename PointMatcher<T>::TransformationParameters TransformationParameters;
	 typedef typename PointMatcher<T>::TransformationCheckers TransformationCheckers;
	 typedef typename PointMatcher<T>::Matrix Matrix;

	inline static const std::string description()
	{
		return "Keep in memory stastistics at each step.";
	}

	struct IterationData {
		TransformationParameters tfParameters;
		DataPoints filteredReference;
		DataPoints reading;
		Matches matches;
		OutlierWeights outlierWeights;
		TransformationCheckers transCheck;
	};

	std::vector<IterationData> iterationsStats;
	std::map<std::string, std::vector<double>> stats;

public:

	virtual void dumpIteration(const size_t iterationNumber, const TransformationParameters& parameters, const DataPoints& filteredReference, const DataPoints& reading, const Matches& matches, const OutlierWeights& outlierWeights, const TransformationCheckers& transformationCheckers);
	virtual void addStat(const std::string& name, double data);

};

} // end of PointMatcherSupport


#endif // __POINTMATCHER_IN_MEMORY_INSPECTORS_H
