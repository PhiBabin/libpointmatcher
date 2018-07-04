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

#ifndef __POINTMATCHER_OUTLIERFILTERS_H
#define __POINTMATCHER_OUTLIERFILTERS_H

#include "PointMatcher.h"

template<typename T>
struct OutlierFiltersImpl
{
	typedef PointMatcherSupport::Parametrizable Parametrizable;
	typedef PointMatcherSupport::Parametrizable P;
	typedef Parametrizable::Parameters Parameters;
	typedef Parametrizable::ParameterDoc ParameterDoc;
	typedef Parametrizable::ParametersDoc ParametersDoc;
	typedef Parametrizable::InvalidParameter InvalidParameter;
	
	typedef typename PointMatcher<T>::DataPoints DataPoints;
	typedef typename PointMatcher<T>::Matches Matches;
	typedef typename PointMatcher<T>::OutlierFilter OutlierFilter;
	typedef typename PointMatcher<T>::InspectorPtr InspectorPtr;
	typedef typename PointMatcher<T>::OutlierWeights OutlierWeights;
	typedef typename PointMatcher<T>::Matrix Matrix;	
	typedef typename PointMatcher<T>::Vector Vector;
	
	struct NullOutlierFilter: public OutlierFilter
	{
		inline static const std::string description()
		{
			return "Does nothing.";
		}
		
		virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
	};

	struct MaxDistOutlierFilter: public OutlierFilter
	{
		inline static const std::string description()
		{
			return "This filter considers as outlier links whose norms are above a fix threshold.";
		}
		inline static const ParametersDoc availableParameters()
		{
			return boost::assign::list_of<ParameterDoc>
				( "maxDist", "threshold distance (Euclidean norm)", "1", "0.0000001", "inf", &P::Comp<T>) 
			;
		}
		
		const T maxDist;
		
		MaxDistOutlierFilter(const Parameters& params = Parameters());
		virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
	};

	//FIXME: is that useful in any case?
	struct MinDistOutlierFilter: public OutlierFilter
	{
		inline static const std::string description()
		{
			return "This filter considers as outlier links whose norms are below a threshold.";
		}
		inline static const ParametersDoc availableParameters()
		{
			return boost::assign::list_of<ParameterDoc>
				( "minDist", "threshold distance (Euclidean norm)", "1", "0.0000001", "inf", &P::Comp<T>) 
			;
		}
		
		const T minDist;
		
		MinDistOutlierFilter(const Parameters& params = Parameters()); 
		virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
	};

	struct MedianDistOutlierFilter: public OutlierFilter 
	{
		inline static const std::string description()
		{
			return "This filter considers as outlier links whose norms are above the median link norms times a factor. Based on \\cite{Diebel2004Median}.";
		}
		inline static const ParametersDoc availableParameters()
		{
			return boost::assign::list_of<ParameterDoc>
				( "factor", "points farther away factor * median will be considered outliers.", "3", "0.0000001", "inf", &P::Comp<T>)
			;
		}
		
		const T factor;
		
		MedianDistOutlierFilter(const Parameters& params = Parameters());
		virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
	};

	struct TrimmedDistOutlierFilter: public OutlierFilter
	{
		inline static const std::string description()
		{
			return "Hard rejection threshold using quantile. This filter considers as inlier a certain percentage of the links with the smallest norms. Based on \\cite{Chetverikov2002Trimmed}.";
		}
		inline static const ParametersDoc availableParameters()
		{
			return boost::assign::list_of<ParameterDoc>
				( "ratio", "percentage to keep", "0.85", "0.0000001", "1.0", &P::Comp<T>)
			;
		}
		
		const T ratio;
		
		TrimmedDistOutlierFilter(const Parameters& params = Parameters());
		virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
	};

	struct VarTrimmedDistOutlierFilter: public OutlierFilter
	{
		inline static const std::string description()
		{
			return "Hard rejection threshold using quantile and variable ratio. Based on \\cite{Phillips2007VarTrimmed}.";
		}
		inline static const ParametersDoc availableParameters()
		{
			return boost::assign::list_of<ParameterDoc>
				( "minRatio", "min ratio", "0.05", "0.0000001", "1", &P::Comp<T>)
				( "maxRatio", "max ratio", "0.99", "0.0000001", "1", &P::Comp<T>)
				( "lambda", "lambda (part of the term that balance the rmsd: 1/ratio^lambda", "0.95" )
			;
		}

		const T minRatio;
		const T maxRatio;
		const T lambda;
		std::vector<T> tunedRatios;
		
		VarTrimmedDistOutlierFilter(const Parameters& params = Parameters());
		virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
		virtual void addStat(InspectorPtr& inspector) const;

	private:
		// return the optimized ratio
		T optimizeInlierRatio(const Matches& matches);
	};
	
	
	// ---------------------------------
	struct SurfaceNormalOutlierFilter: public OutlierFilter
	{
		inline static const std::string description()
		{
			return "Hard rejection threshold using the angle between the surface normal vector of the reading and the reference. If normal vectors or not in the descriptor for both of the point clouds, does nothing.";
		}
		inline static const ParametersDoc availableParameters()
		{
			return boost::assign::list_of<ParameterDoc>
				( "maxAngle", "Maximum authorised angle between the 2 surface normals (in radian)", "1.57", "0.0", "3.1416", &P::Comp<T>)
			;
		}

		const T eps;
		bool warningPrinted;
		
		SurfaceNormalOutlierFilter(const Parameters& params = Parameters());
		virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
	};

	struct GenericDescriptorOutlierFilter: public OutlierFilter
	{
		inline static const std::string description()
		{
			return "This filter weights matched points based on a 1D descriptor of either a single point cloud (either the reference or the reading). The descriptor values must be larger than zero.";
		}
		inline static const ParametersDoc availableParameters()
		{
			return boost::assign::list_of<ParameterDoc>
				( "source", "Point cloud from which the descriptor will be used: reference or reading", "reference")
				( "descName", "Descriptor name used to weight paired points", "none")
				( "useSoftThreshold", "If set to 1 (true), uses the value of the descriptor as a weight. If set to 0 (false), uses the parameter 'threshold' to set binary weights.", "0", "0", "1", P::Comp<bool>)
				( "useLargerThan", "If set to 1 (true), values over the 'threshold' will have a weight of one.  If set to 0 (false), values under the 'threshold' will have a weight of one. All other values will have a weight of zero.", "1", "0", "1", P::Comp<bool>)
				( "threshold", "Value used to determine the binary weights", "0.1", "0.0000001", "inf", &P::Comp<T>)
				;
		}
		
		const std::string source;
		const std::string descName;
		const bool useSoftThreshold;
		const bool useLargerThan;
		const T threshold;
		
		GenericDescriptorOutlierFilter(const Parameters& params = Parameters());
		virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
	};

    struct RobustOutlierFilter: public OutlierFilter
    {

		inline static const std::string description()
		{
            return "Robust weight function part of the M-Estimator familly. 7 robust functions to choose from (Cauchy, Welsch, Switchable Constraint, Geman-McClure, Tukey, Huber and L1). The Welsch weight uses an exponential decay reducing the influence of matched point farther away \\cite{RobustWeightFunctions}. More explicitly, the function is w = exp[- (matched distance)^2/tuning^2].";
		}
		inline static const ParametersDoc availableParameters()
		{
            return boost::assign::list_of<ParameterDoc>
                ( "robustFct", "Type of robust function used. Available fct: 'cauchy', 'welsch', 'sc'(aka Switchable-Constraint), 'gm' (aka Geman-McClure), 'tukey', 'huber' and 'L1'. (Default: cauchy)", "cauchy")
                ( "scale", "Tuning parameter used to limit the influence of outliers. It could be interpreted as a standard deviation. The unit of this parameter is the same as the distance used, typically meters.", "0.2", "0.0000001", "inf", &P::Comp<T>)
				        ( "approximation", "If the matched distance is larger than this threshold, its weight will be forced to zero. This can save computation as zero values are not minimized. If set to inf (default value), no approximation is done. The unit of this parameter is the same as the distance used, typically meters.", "inf", "0.0", "inf", &P::Comp<T>)
                ( "useMadForScale", "Instead of using the tuning parameter, the Median of Absolute Deviations(MAD) is used to tune the influence of outliers.", "0", "0", "1", &P::Comp<bool>)
                ( "nbMadIteration", "For how many iteration MAD is recalculate. This parameter is only use if useMadForScale==1. If nbMadIteration==0 then it's recalculate at every iteration", "0", "0", "100", &P::Comp<int>)
                ( "useBergstromScale", "Use an iterative exponentialy decreasive value for the tuning", "0", "0", "1", &P::Comp<bool>)
                    ;
		}


        virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
      RobustOutlierFilter(const std::string& className, const ParametersDoc paramsDoc, const Parameters& params);
      RobustOutlierFilter(const Parameters& params = Parameters());
		    virtual void addStat(InspectorPtr& inspector) const;
    protected:
			enum RobustFctId {
					Cauchy=0,
					Welsch=1,
					SwitchableConstraint=2,
					GM=3,
					Tukey=4,
					Huber=5,
					L1=6,
          Student=7,
					Lp=8
			};
			typedef std::map<std::string, RobustFctId> RobustFctMap;
      static RobustFctMap robustFcts;
      const std::string robustFctName;
      T tuning;
      const T squaredApproximation;
      const bool useMad;
      const int nbMadIteration;
      const bool useBergstromScale;
      int robustFctId;
      int iteration;
      T scale;


      virtual void resolveEstimatorName();
      virtual OutlierWeights robustFiltering(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
	};

    struct RobustTrimmedOutlierFilter: public RobustOutlierFilter
    {

      inline static const std::string description()
      {
        return "Like RobustOutlierFilter, but use a trimmed distance filter as a final filtering step.";
      }
      inline static const ParametersDoc availableParameters()
      {
        return boost::assign::list_of<ParameterDoc>
                ( "robustFct", "Type of robust function used. Available fct: 'cauchy', 'welsch', 'sc'(aka Switchable-Constraint), 'gm' (aka Geman-McClure), 'tukey', 'huber' and 'L1'. (Default: cauchy)", "cauchy")
                ( "ratio", "Percentile of the distance to keep, this is an estimation of the overlap between the scans", "0.2", "0.0000001", "inf", &P::Comp<T>)
                ( "scale", "Tuning parameter used to limit the influence of outliers. It could be interpreted as a standard deviation. The unit of this parameter is the same as the distance used, typically meters.", "0.2", "0.0000001", "inf", &P::Comp<T>)
                ( "approximation", "If the matched distance is larger than this threshold, its weight will be forced to zero. This can save computation as zero values are not minimized. If set to inf (default value), no approximation is done. The unit of this parameter is the same as the distance used, typically meters.", "inf", "0.0", "inf", &P::Comp<T>)
                ( "useMadForScale", "Instead of using the tuning parameter, the Median of Absolute Deviations(MAD) is used to tune the influence of outliers.", "0", "0", "1", &P::Comp<bool>)
                ( "nbMadIteration", "For how many iteration MAD is recalculate. This parameter is only use if useMadForScale==1. If nbMadIteration==0 then it's recalculate at every iteration", "0", "0", "100", &P::Comp<int>)
                ( "useBergstromScale", "Use an iterative exponentialy decreasive value for the tuning", "0", "0", "1", &P::Comp<bool>)
                ;
      }

      const T ratio;

      virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
      RobustTrimmedOutlierFilter(const Parameters& params = Parameters());
    };

    struct IterativeCauchyOutlierFilter: public OutlierFilter
    {
        inline static const std::string description()
        {
          return "Bosse's implementation of Cauchy, it multiply the tuning constant by the iteration number";
        }
        inline static const ParametersDoc availableParameters()
        {
          return boost::assign::list_of<ParameterDoc>
                  ( "scale", "Tuning parameter used to limit the influence of outliers. It could be interpreted as a standard deviation. The unit of this parameter is the same as the distance used, typically meters.", "0.2", "0.0000001", "inf", &P::Comp<T>)
                  ;
        }

        const T scale;
        int iteration;

        IterativeCauchyOutlierFilter(const Parameters& params = Parameters());
        virtual OutlierWeights compute(const DataPoints& filteredReading, const DataPoints& filteredReference, const Matches& input);
    };

}; // OutlierFiltersImpl

#endif // __POINTMATCHER_OUTLIERFILTERS_H
