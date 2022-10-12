/*
 * StressTuple.h
 *
 *  Created on: Dec 15, 2019
 *      Author: Nikhil
 */

#ifndef INCLUDE_STRESSTUPLE_H_
#define INCLUDE_STRESSTUPLE_H_

// Recursively build stress of a type sType= Piola/Cauchy
template<std::size_t I=0, typename ...TStress>
inline typename std::enable_if<I == sizeof...(TStress), void>::type
 recursiveBuildStress(const double& fij,
					  const Vector3d& ra,
					  const Vector3d& rA,
					  const Vector3d& rb,
					  const Vector3d& rB,
					  const Vector3d& rab,
					  const Vector3d& rAB,
					  const int& i_gridPoint,
					  const int& i_stress,
					  std::tuple<TStress&...> t)
{
	if (sizeof...(TStress)!=0)
		assert(0);
}
template<std::size_t I=0, StressType stressType, typename ...BF>
inline typename std::enable_if<I < sizeof...(BF), void>::type
 recursiveBuildStress(const double& fij,
					  const Vector3d& ra,
					  const Vector3d& rA,
					  const Vector3d& rb,
					  const Vector3d& rB,
					  const Vector3d& rab,
					  const Vector3d& rAB,
					  const int& i_gridPoint,
					  const int& i_stress,
					  std::tuple<Stress<BF,stressType>&...> t)
{
	if (I == i_stress && stressType == Cauchy)
	{
		assert(rab.squaredNorm()>epsilon);
		std::get<I>(t).field[i_gridPoint]= std::get<I>(t).field[i_gridPoint] +
				std::get<I>(t).method.bondFunction(ra,rb)*fij*rab.transpose()*rab/rab.norm();
	}
	else if (I == i_stress && stressType == Piola)
	{
		assert(rab.squaredNorm()>epsilon);
		std::get<I>(t).field[i_gridPoint]= std::get<I>(t).field[i_gridPoint] +
				std::get<I>(t).method.bondFunction(rA,rB)*fij*rAB.transpose()*rab/rab.norm();
	}
	else
		recursiveBuildStress<I+1>(fij,ra,rA,rb,rB,rab,rAB,i_gridPoint,i_stress,t);
}

// Recursively nullify stress of a type sType= Piola/Cauchy
template<std::size_t I=0, typename ...TStress>
inline typename std::enable_if<I == sizeof...(TStress), void>::type
recursiveNullifyStress(std::tuple<TStress&...> t)
{
}
template<std::size_t I=0, StressType stressType, typename ...BF>
inline typename std::enable_if<I < sizeof...(BF), void>::type
recursiveNullifyStress(std::tuple<Stress<BF,stressType>&...> t)
{
    std::fill(std::get<I>(t).field.begin(), std::get<I>(t).field.end(),Matrix3d::Zero() );
    recursiveNullifyStress<I+1>(t);
}

//////////  Get the maximum averaging domain size across all bond functions of stresses of type stressType ////////////
double averagingDomainSize_max(const std::tuple<> t)
{
	return 0;
}
template<std::size_t I=0, typename ...TStress>
inline typename std::enable_if<I == sizeof...(TStress)-1, double>::type
	averagingDomainSize_max(const std::tuple<TStress&...> t)
{
	return std::get<I>(t).method.getAveragingDomainSize();
}
template<std::size_t I=0, typename ...TStress>
inline typename std::enable_if<I < sizeof...(TStress)-1, double>::type
 averagingDomainSize_max(const std::tuple<TStress&...> t)
{
	return std::max(std::get<I>(t).method.getAveragingDomainSize(),averagingDomainSize_max<I+1>(t));
}

///////////////////////////      Get the maximum averaging domain sizes for each grid     //////////////////////////////
///////////////////////////	 Returns a mapping between grid pointers and maximum domain size   /////////////////////////
std::map<Grid<Reference>*,double> recursiveGridMaxAveragingDomainSizeMap(const std::tuple<>& t)
{
std::map<Grid<Reference>*,double> map;
return map;
}
template<std::size_t I=0, StressType stressType,
		 typename TGrid= typename std::conditional<stressType == Piola,Grid<Reference>,Grid<Current>>::type,
		 typename ...BF>
inline typename std::enable_if<I == sizeof...(BF), std::map<TGrid*,double>>::type
 recursiveGridMaxAveragingDomainSizeMap(const std::tuple<Stress<BF,stressType>&...> t)
{
	std::map<TGrid*,double> map;
	return map;
}
template<std::size_t I=0, StressType stressType,
		 typename TGrid= typename std::conditional<stressType == Piola,Grid<Reference>,Grid<Current>>::type,
		 typename ...BF>
inline typename std::enable_if<I < sizeof...(BF), std::map<TGrid*,double>>::type
 recursiveGridMaxAveragingDomainSizeMap(const std::tuple<Stress<BF,stressType>&...> t)
{
	std::map<TGrid*,double> map;
	const std::map<TGrid*,double>& rmap= recursiveGridMaxAveragingDomainSizeMap<I+1>(t);
	map[std::get<I>(t).pgrid]= std::get<I>(t).method.getAveragingDomainSize();

	// Loop over rmap
	for (const auto& pair : rmap)
		if (std::get<I>(t).pgrid == pair.first)
		{
			map[pair.first]= std::max(map.at(std::get<I>(t).pgrid),pair.second);
		}
		else
		{
			auto result= map.insert(pair);
			assert(result.second==true);
		}
	return map;
}

///// Returns a vector of pairs, where pair = (grid pointer, averaging domain size) for a given type= Reference/Current /////////////

inline std::vector<std::pair<Grid<Reference>*,double>> getTGridDomainSizePairs(const std::tuple<> emptyTuple)
{
	std::vector<std::pair<Grid<Reference>*,double>> emptyVectorPair;
	return emptyVectorPair;
}
template<std::size_t I=0, StressType stressType,
		 typename TGrid= typename std::conditional<stressType == Piola,Grid<Reference>,Grid<Current>>::type,
		 typename ...BF>
inline typename std::enable_if<I == sizeof...(BF)-1, std::vector<std::pair<TGrid*,double>>>::type
 getTGridDomainSizePairs(const std::tuple<Stress<BF,stressType>&...> t)
{
	std::vector<std::pair<TGrid*,double>> vectorPair;
	vectorPair.push_back({std::get<I>(t).pgrid,std::get<I>(t).method.getAveragingDomainSize()});
	return vectorPair;
}
template<std::size_t I=0, StressType stressType,
		 typename TGrid= typename std::conditional<stressType == Piola,Grid<Reference>,Grid<Current>>::type,
 		 typename ...BF>
inline typename std::enable_if< I < sizeof...(BF)-1, std::vector<std::pair<TGrid*,double>> >::type
 getTGridDomainSizePairs(const std::tuple<Stress<BF,stressType>&...> t)
{
	std::vector<std::pair<TGrid*,double>> vectorPair;
	vectorPair.push_back({std::get<I>(t).pgrid,std::get<I>(t).method.getAveragingDomainSize()});
	std::vector<std::pair<TGrid*,double>> next= getTGridDomainSizePairs<I+1>(t);
	vectorPair.insert(vectorPair.end(),next.begin(),next.end());
	return vectorPair;
}

/////////////// Returns a vector of all grid pointers for type= Reference/Current /////////////////
std::vector<GridBase*> getBaseGridList(const std::tuple<> t)
{
	std::vector<GridBase*> pgridBaseVector;
	return pgridBaseVector;
}
template<std::size_t I=0, typename ...TStress>
inline typename std::enable_if<I == sizeof...(TStress)-1, std::vector<GridBase*>>::type
 getBaseGridList(const std::tuple<TStress&...> t)
{
	std::vector<GridBase*> pgridBaseVector;
	pgridBaseVector.push_back(std::get<I>(t).pgrid);
	return pgridBaseVector;
}
template<std::size_t I=0,typename ...TStress>
inline typename std::enable_if<I < sizeof...(TStress)-1, std::vector<GridBase*>>::type
 getBaseGridList(const std::tuple<TStress&...> t)
{
	std::vector<GridBase*> pgridBaseVector;

	pgridBaseVector.push_back(std::get<I>(t).pgrid);
	std::vector<GridBase*> next= getBaseGridList<I+1>(t);
	pgridBaseVector.insert(pgridBaseVector.end(),next.begin(),next.end());
	return pgridBaseVector;

}

///////////////////////////      Fold the grid points in all the grids of type= Reference/Current     //////////////////////////////
template<std::size_t I=0, typename ...TStress>
inline typename std::enable_if<I == sizeof...(TStress), void>::type
 recursiveFold(const Vector3d& origin, const Vector3d& step, const Vector3i& pbc, const std::tuple<TStress&...> t)
{ }
template<std::size_t I=0, typename ...TStress>
inline typename std::enable_if<I < sizeof...(TStress), void>::type
 recursiveFold(const Vector3d& origin,
						  const Vector3d& step,
						  const Vector3i& pbc,
						  std::tuple<TStress&...> t)
{
	BoxPoints boxPoints(origin,step, std::get<I>(t).pgrid->coordinates);
	std::cout << "Folding grid points in grid: " << std::get<I>(t).pgrid << " if necessary...."<<"\n";
	boxPoints.fold(pbc);
	std::cout << std::endl;
	recursiveFold<I+1>(origin,step,pbc,t);
}

// Write grid data for each requested stress field
template<std::size_t I=0, typename ...TStress>
inline typename std::enable_if<I == sizeof...(TStress), void>::type
 recursiveWriteStressAndGrid(std::tuple<TStress&...> t)
{ }

template<std::size_t I=0, typename ...TStress>
inline typename std::enable_if<I < sizeof...(TStress), void>::type
 recursiveWriteStressAndGrid(std::tuple<TStress&...> t)
{
	std::cout << "Writing stress and grid of " << std::get<I>(t).name << std::endl;
	std::get<I>(t).write();
	std::get<I>(t).pgrid->write(std::get<I>(t).name);
	recursiveWriteStressAndGrid<I+1>(t);
}

#endif

