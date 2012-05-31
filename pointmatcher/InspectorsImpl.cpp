// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Copyright (c) 2010--2011,
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

#include "InspectorsImpl.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

template<typename T>
InspectorsImpl<T>::PerformanceInspector::PerformanceInspector(const std::string& className, const ParametersDoc paramsDoc, const Parameters& params):
	Inspector(className,paramsDoc,params)
{
}

template<typename T>
InspectorsImpl<T>::PerformanceInspector::PerformanceInspector(const Parameters& params):
	Inspector("PerformanceInspector", PerformanceInspector::availableParameters(), params)
{}

template<typename T>
void InspectorsImpl<T>::PerformanceInspector::addStat(const std::string& name, double data)
{
	HistogramMap::iterator it(stats.find(name));
	if (it == stats.end())
		it = stats.insert(
			HistogramMap::value_type(name, 
				Histogram(16, name, Parametrizable::get<string>("baseFileName"), Parametrizable::get<bool>("dumpPerfOnExit"))
			)
		).first;
	it->second.push_back(data);
}

template<typename T>
void InspectorsImpl<T>::PerformanceInspector::dumpStats(std::ostream& stream)
{
	// Note: this dump format will most probably change in the future
	for (auto it(stats.begin()); it != stats.end(); ++it)
	{
		it->second.dumpStats(stream);
		auto jt(it);
		++jt;
		if (jt != stats.end())
			stream << ", ";
	}
}

template<typename T>
void InspectorsImpl<T>::PerformanceInspector::dumpStatsHeader(std::ostream& stream)
{
	for (auto it(stats.begin()); it != stats.end(); ++it)
	{
		it->second.dumpStatsHeader(stream);
		auto jt(it);
		++jt;
		if (jt != stats.end())
			stream << ", ";
	}
}

template struct InspectorsImpl<float>::PerformanceInspector;
template struct InspectorsImpl<double>::PerformanceInspector;

/*
	// This code creates a header according to information from the datapoints
	std::ofstream ofs(fileName.c_str());
	size_t labelIndex(0);
	size_t labelComponent(0);
	const typename DataPoints::Labels labels(data.featureLabels);
	for (int i = 0; i < features.rows(); ++i)
	{
		// display label for row
		if (labelIndex <= labels.size())
		{
			const typename DataPoints::Label label(labels[labelIndex]);
			ofs << label.text;
			if (label.span > 1)
				ofs << "_" << labelComponent;
			++labelComponent;
			if (labelComponent >= label.span)
				++labelIndex;
		}
		else
		{
			ofs << "?";
		}
		ofs << " ";
	}*/
template<typename T>
InspectorsImpl<T>::AbstractVTKInspector::AbstractVTKInspector(const std::string& className, const ParametersDoc paramsDoc, const Parameters& params):
	PerformanceInspector(className,paramsDoc,params),
	streamIter(0)
{
}

	
template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::dumpDataPoints(const DataPoints& data, std::ostream& stream)
{
	const typename DataPoints::Features& features(data.features);
	//const typename DataPoints::Descriptors& descriptors(data.descriptors);
	
	stream << "# vtk DataFile Version 3.0\n";
	stream << "File created by libpointmatcher\n";
	stream << "ASCII\n";
	stream << "DATASET POLYDATA\n";
	stream << "POINTS " << features.cols() << " float\n";
	if(features.rows() == 4)
	{
		stream << features.topLeftCorner(3, features.cols()).transpose() << "\n";
	}
	else
	{
		stream << features.transpose() << "\n";
	}
	
	stream << "VERTICES "  << features.cols() << " "<< features.cols() * 2 << "\n";
	for (int i = 0; i < features.cols(); ++i)
		stream << "1 " << i << "\n";
	
	stream << "POINT_DATA " << features.cols() << "\n";
	buildScalarStream(stream, "densities", data);
	buildScalarStream(stream, "obstacles", data);
	buildScalarStream(stream, "inclination", data);
	buildScalarStream(stream, "groupId", data);
	buildScalarStream(stream, "simpleSensorNoise", data);
	buildNormalStream(stream, "normals", data);
	buildVectorStream(stream, "eigValues", data);
	buildTensorStream(stream, "eigVectors", data);
	//TODO: add support for colors

}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::dumpMeshNodes(const DataPoints& data, std::ostream& stream)
{
	//const typename DataPoints::Features& features(data.features);
	const typename DataPoints::Descriptors& descriptors(data.descriptors.transpose());
	
	assert(descriptors.cols() >= 15);

	stream << "# vtk DataFile Version 3.0\n";
	stream << "Triangle mesh\n";
	stream << "ASCII\n";
	stream << "DATASET POLYDATA\n";

	stream << "POINTS " << descriptors.rows() * 3 << " float\n"; // not optimal: points and edges are drawn several times!
	for (int i = 0; i < descriptors.rows(); i++)
	{
		// TODO: use getDescriptorByName(nameTag) to access blocks 
		stream << descriptors.block(i, 3, 1, 3) << "\n";
		stream << descriptors.block(i, 6, 1, 3) << "\n";
		stream << descriptors.block(i, 9, 1, 3) << "\n";
	}

	//TODO: add centroids...
	/*
	stream << features.transpose() << "\n";
	*/

	stream << "POLYGONS " << descriptors.rows() << " " << descriptors.rows() * 4 << "\n";
	for (int i = 0; i < descriptors.rows(); i++)
	{
		stream << "3 " << (i*3) << " " << (i*3 + 1) << " " << (i*3 + 2) << "\n";
	}
	
	stream << "CELL_DATA " << descriptors.rows() << "\n";

	stream << "NORMALS triangle_normals float\n";
	stream << descriptors.block(0, 0, descriptors.rows(), 3) << "\n";

	// TODO: use descriptor labels, particularly to represent normals
	//stream << "POINT_DATA " << descriptors.rows() << "\n";
	//buildScalarStream(stream, "densities", data);
	//buildNormalStream(stream, "normals", data);
	//buildVectorStream(stream, "eigValues", data);
	//buildTensorStream(stream, "eigVectors", data);
}

// FIXME:rethink how we dump stuff (accumulate in a correctly-referenced table, and then dump?) and unify with previous
template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::dumpDataLinks(
	const DataPoints& ref, 
	const DataPoints& reading, 
	const Matches& matches, 
	const OutlierWeights& featureOutlierWeights, 
	std::ostream& stream)
{

	const typename DataPoints::Features& refFeatures(ref.features);
	const int refPtCount(refFeatures.cols());
	//const int featDim(refFeatures.rows());
	const typename DataPoints::Features& readingFeatures(reading.features);
	const int readingPtCount(readingFeatures.cols());
	const int totalPtCount(refPtCount+readingPtCount);
	
	stream << "# vtk DataFile Version 3.0\n";
	stream << "comment\n";
	stream << "ASCII\n";
	stream << "DATASET POLYDATA\n";
	
	stream << "POINTS " << totalPtCount << " float\n";
	if(refFeatures.rows() == 4)
	{
		// reference pt
		stream << refFeatures.topLeftCorner(3, refFeatures.cols()).transpose() << "\n";
		// reading pt
		stream << readingFeatures.topLeftCorner(3, readingFeatures.cols()).transpose() << "\n";
	}
	else
	{
		// reference pt
		stream << refFeatures.transpose() << "\n";
		// reading pt
		stream << readingFeatures.transpose() << "\n";
	}
	const int knn = matches.ids.rows();
	
	stream << "LINES " << readingPtCount*knn << " "  << readingPtCount*knn * 3 << "\n";
	//int j = 0;
	for(int k = 0; k < knn; k++) // knn
	{
		for (int i = 0; i < readingPtCount; ++i)
		{
			stream << "2 " << refPtCount + i << " " << matches.ids(k, i) << "\n";
		}
	}

	stream << "CELL_DATA " << readingPtCount*knn << "\n";
	stream << "SCALARS outlier float 1\n";
	stream << "LOOKUP_TABLE default\n";
	//stream << "LOOKUP_TABLE alphaOutlier\n";
	for(int k = 0; k < knn; k++) // knn
	{
		for (int i = 0; i < readingPtCount; ++i) //nb pts
		{
			stream << featureOutlierWeights(k, i) << "\n";
		}
	}

	//stream << "LOOKUP_TABLE alphaOutlier 2\n";
	//stream << "1 0 0 0.5\n";
	//stream << "0 1 0 1\n";

}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::dumpDataPoints(const DataPoints& filteredReference, const std::string& name)
{
	ostream* stream(openStream(name));
	dumpDataPoints(filteredReference, *stream);
	closeStream(stream);
}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::dumpMeshNodes(const DataPoints& filteredReference, const std::string& name)
{
	ostream* stream(openStream(name));
	dumpMeshNodes(filteredReference, *stream);
	closeStream(stream);
}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::dumpIteration(
	const size_t iterationCount,
	const TransformationParameters& parameters,
	const DataPoints& filteredReference,
	const DataPoints& reading,
	const Matches& matches,
	const OutlierWeights& outlierWeights, 
	const TransformationCheckers& transCheck)
{
	ostream* streamLinks(openStream("link", iterationCount));
	dumpDataLinks(filteredReference, reading, matches, outlierWeights, *streamLinks);
	closeStream(streamLinks);

	ostream* streamRead(openStream("reading", iterationCount));
	dumpDataPoints(reading, *streamRead);
	closeStream(streamRead);
	
	ostream* streamRef(openStream("reference", iterationCount));
	dumpDataPoints(filteredReference, *streamRef);
	closeStream(streamRef);
	
	// streamIter must be define by children
	assert(streamIter);

	if(iterationCount == 0)
	{
		//Build header
		for(unsigned int j = 0; j < transCheck.size(); j++)
		{
			for(unsigned int i=0; i < transCheck[j]->getValueNames().size(); i++)
			{
				if (!(j == 0 && i == 0))
					*streamIter << ", ";
				*streamIter << transCheck[j]->getValueNames()[i] << ", "; 
				*streamIter << transCheck[j]->getLimitNames()[i]; 
			}
		}
		
		*streamIter << "\n";
	}

	for(unsigned int j = 0; j < transCheck.size(); j++)
	{
		for(unsigned int i=0; i < transCheck[j]->getValueNames().size(); i++)
		{
		
			if (!(j == 0 && i == 0))
				*streamIter << ", ";
			*streamIter << transCheck[j]->getValues()(i) << ", ";
			*streamIter << transCheck[j]->getLimits()(i); 
		}
	}

	*streamIter << "\n";
}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::buildGenericAttributeStream(std::ostream& stream, const std::string& attribute, const std::string& nameTag, const DataPoints& cloud, const int forcedDim)
{
	const auto desc(cloud.getDescriptorViewByName(nameTag));	

	assert(desc.rows() <= forcedDim);

	if(desc.rows() != 0 && desc.rows() != 0)
	{
		stream << attribute << " " << nameTag << " float\n";
		if(attribute.compare("SCALARS") == 0)
			stream << "LOOKUP_TABLE default\n";

		stream << padWithZeros(desc, forcedDim, desc.cols()).transpose();
		stream << "\n";
	}
}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::buildScalarStream(std::ostream& stream,
	const std::string& name,
	const DataPoints& cloud)
{
	buildGenericAttributeStream(stream, "SCALARS", name, cloud, 1);
}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::buildNormalStream(std::ostream& stream,
	const std::string& name,
	const DataPoints& cloud)
{
	buildGenericAttributeStream(stream, "NORMALS", name, cloud, 3);
}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::buildVectorStream(std::ostream& stream,
	const std::string& name,
	const DataPoints& cloud)
{
	buildGenericAttributeStream(stream, "VECTORS", name, cloud, 3);
}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::buildTensorStream(std::ostream& stream,
	const std::string& name,
	const DataPoints& cloud)
{
	buildGenericAttributeStream(stream, "TENSORS", name, cloud, 9);
}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::buildScalarStream(std::ostream& stream,
	const std::string& name,
	const DataPoints& ref, 
	const DataPoints& reading)
{
			
	const Matrix descRef(ref.getDescriptorByName(name));	
	const Matrix descRead(reading.getDescriptorByName(name));

	if(descRef.rows() != 0 && descRead.rows() != 0)
	{
		stream << "SCALARS " << name << " float\n";
		stream << "LOOKUP_TABLE default\n";
		
		stream << padWithZeros(
				descRef, 1, ref.descriptors.cols()).transpose();
		stream << "\n";
		stream << padWithZeros(
				descRead, 1, reading.descriptors.cols()).transpose();
		stream << "\n";
	}
}


template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::buildNormalStream(std::ostream& stream,
	const std::string& name,
	const DataPoints& ref, 
	const DataPoints& reading)
{
			
	const Matrix descRef(ref.getDescriptorByName(name));	
	const Matrix descRead(reading.getDescriptorByName(name));

	if(descRef.rows() != 0 && descRead.rows() != 0)
	{
		stream << "NORMALS " << name << " float\n";

		stream << padWithZeros(
				descRef, 3, ref.descriptors.cols()).transpose();
		stream << "\n";
		stream << padWithZeros(
				descRead, 3, reading.descriptors.cols()).transpose();
		stream << "\n";
	}
}


template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::buildVectorStream(std::ostream& stream,
	const std::string& name,
	const DataPoints& ref, 
	const DataPoints& reading)
{
			
	const Matrix descRef(ref.getDescriptorByName(name));	
	const Matrix descRead(reading.getDescriptorByName(name));

	if(descRef.rows() != 0 && descRead.rows() != 0)
	{
		stream << "VECTORS " << name << " float\n";

		stream << padWithZeros(
				descRef, 3, ref.descriptors.cols()).transpose();
		stream << "\n";
		stream << padWithZeros(
				descRead, 3, reading.descriptors.cols()).transpose();
		stream << "\n";
	}
}


template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::buildTensorStream(std::ostream& stream,
	const std::string& name,
	const DataPoints& ref, 
	const DataPoints& reading)
{
			
	const Matrix descRef(ref.getDescriptorByName(name));	
	const Matrix descRead(reading.getDescriptorByName(name));

	if(descRef.rows() != 0 && descRead.rows() != 0)
	{
		stream << "TENSORS " << name << " float\n";

		stream << padWithZeros(
				descRef, 9, ref.descriptors.cols()).transpose();
		stream << "\n";
		stream << padWithZeros(
				descRead, 9, reading.descriptors.cols()).transpose();
		stream << "\n";
	}
}

template<typename T>
typename PointMatcher<T>::Matrix InspectorsImpl<T>::AbstractVTKInspector::padWithZeros(
	const Matrix m,
	const int expectedRow,
	const int expectedCols)
{
	assert(m.cols() <= expectedCols || m.rows() <= expectedRow);
	if(m.cols() == expectedCols && m.rows() == expectedRow)
	{
		return m;
	}
	else
	{
		Matrix tmp = Matrix::Zero(expectedRow, expectedCols); 
		tmp.topLeftCorner(m.rows(), m.cols()) = m;
		return tmp;
	}

}

template<typename T>
void InspectorsImpl<T>::AbstractVTKInspector::finish(const size_t iterationCount)
{
}


//-----------------------------------
// VTK File inspector

template<typename T>
InspectorsImpl<T>::VTKFileInspector::VTKFileInspector(const Parameters& params):
	AbstractVTKInspector("VTKFileInspector", VTKFileInspector::availableParameters(), params),
	baseFileName(Parametrizable::get<string>("baseFileName"))
{
}

template<typename T>
void InspectorsImpl<T>::VTKFileInspector::init()
{
	ostringstream oss;
	oss << baseFileName << "-iterationInfo.csv";

	this->streamIter = new ofstream(oss.str().c_str());
	if (this->streamIter->fail())
		throw std::runtime_error("Couldn't open the file \"" + oss.str() + "\". Check if repository exist.");
	
}

template<typename T>
void InspectorsImpl<T>::VTKFileInspector::finish(const size_t iterationCount)
{
	closeStream(this->streamIter);
}

template<typename T>
std::ostream* InspectorsImpl<T>::VTKFileInspector::openStream(const std::string& role)
{
	string filteredStr = role;
	if(role.substr(role.size()-4,4) == ".vtk")
		filteredStr = role.substr(0, role.size()-4);

	ostringstream oss;
	if(baseFileName != "")
		oss << baseFileName << "-" << filteredStr << ".vtk";
	else
		oss << filteredStr << ".vtk";

	ofstream* file = new ofstream(oss.str().c_str());
	if (file->fail())
		throw std::runtime_error("Couldn't open the file \"" + oss.str() + "\". Check if repository exist.");
	return file;
}

template<typename T>
std::ostream* InspectorsImpl<T>::VTKFileInspector::openStream(const std::string& role, const size_t iterationCount)
{
	ostringstream oss;
	oss << baseFileName << "-" << role << "-" << iterationCount << ".vtk";
	ofstream* file = new ofstream(oss.str().c_str());
	if (file->fail())
		throw std::runtime_error("Couldn't open the file \"" + oss.str() + "\". Check if repository exist.");
	return file;
}

template<typename T>
void InspectorsImpl<T>::VTKFileInspector::closeStream(std::ostream* stream)
{
	delete stream;
}



template struct InspectorsImpl<float>::VTKFileInspector;
template struct InspectorsImpl<double>::VTKFileInspector;
