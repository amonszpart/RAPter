#ifndef __GfxTL_NULLSTRATEGY_HPP__
#define __GfxTL_NULLSTRATEGY_HPP__

namespace GfxTL
{

	template< class Base >
	void NullStrategy::StrategyBase< Base >::StartPointInfoCache()
	{}

	template< class Base >
	void NullStrategy::StrategyBase< Base >::StopPointInfoCache()
	{}

	template< class Base >
	bool NullStrategy::StrategyBase< Base >::ShouldSubdivide(
		const CellType &cell) const
	{
		return false;
	}

	template< class Base >
	void NullStrategy::StrategyBase< Base >::InitCellData(CellType *)
	{}

	template< class Base >
	void NullStrategy::StrategyBase< Base >::
		InitLevelDependentCellData(CellType *)
	{}

};

#endif // __GfxTL_NULLSTRATEGY_HPP__
