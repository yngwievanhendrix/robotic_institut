
/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_VectorND.h
 *  \brief  Definition of class \c GW_VectorND
 *  \author Gabriel Peyré
 *  \date   5-31-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _GW_VECTORND_H_
#define _GW_VECTORND_H_


//#error toto
#include "GW_MathsConfig.h"
#include "GW_VectorStatic.h"

#include "tnt/tnt.h"
#include "tnt/jama_eig.h"
#include "tnt/jama_qr.h"
#include "tnt/jama_cholesky.h"
#include "tnt/jama_svd.h"
#include "tnt/jama_lu.h"


namespace GW {

typedef TNT::Array1D<GW_Float> GW_TNTArray1D;

/*------------------------------------------------------------------------------*/
/** 
 *  \class  GW_VectorND
 *  \brief  A vector of size N
 *  \author Gabriel Peyré
 *  \date   5-31-2003
 *
 *  Use \b TNT lib
 */ 
/*------------------------------------------------------------------------------*/

class GW_VectorND:		public GW_TNTArray1D 
{

public:

    /*------------------------------------------------------------------------------*/
    /** \name Constructor and destructor */
    /*------------------------------------------------------------------------------*/
    //@{
	GW_VectorND()							:GW_TNTArray1D(){}
	explicit GW_VectorND(GW_U32 n)			:GW_TNTArray1D((int) n){}
	GW_VectorND(GW_U32 n,  GW_Float *a)		:GW_TNTArray1D((int) n,a){}
	GW_VectorND(GW_U32 n, const GW_Float &a):GW_TNTArray1D((int) n,a){}
	GW_VectorND(const GW_VectorND &A)		:GW_TNTArray1D(A){}
    //@}

	GW_U32 GetDim() const	
	{ 
		return this->dim(); 
	}

	void Reset( GW_U32 dim )
	{
		((GW_TNTArray1D&) (*this)) = GW_TNTArray1D((int) dim);
	}

	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  operator
	*
	*  \param  f Rigth hand statement
	*  \return the vector scaled.
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        GW_VectorND operator*(const GW_Float& f) const
	{
		GW_VectorND Res( this->GetDim() );
                 Multiply( f, *this, Res );
		return Res;
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  Multiply
	*
	*  \param  f left hand statement.
	*  \param  v right hand statement.
	*  \param  r vector v scaled by f.
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        static void Multiply(const GW_Float f, const GW_VectorND& v, GW_VectorND& r)
	{
		GW_ASSERT( v.GetDim() == r.GetDim() );

		for( GW_U32 i=0; i<r.GetDim(); ++i )
			r.SetData( i, f*v.GetData(i) );
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  *=
	*
	*  \param  f scale factOML
	*  \author Gabriel Peyré 2001-09-29
	*
	*	Auto-scale the vector.
	*/ 
	/*------------------------------------------------------------------------------*/
        void operator*=(const GW_Float & f)
	{	
                 Multiply( f, *this, *this );
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  +
	*
	*  \param  v right hand statement
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        GW_VectorND operator+(const GW_VectorND& v) const
	{
		GW_VectorND Res( this->GetDim() );
                 Add( *this, v, Res );
		return Res;
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  Add
	*
	*  \param  a left hand statement
	*  \param  b right hand statement
	*  \param  r result
	*  \author Gabriel Peyré 2001-09-29
	*
	*	Add the two vectOMLs.
	*/ 
	/*------------------------------------------------------------------------------*/
        static void Add(const GW_VectorND& a, const GW_VectorND& b, GW_VectorND& r)
	{
		GW_ASSERT( a.GetDim() == r.GetDim() );
		GW_ASSERT( b.GetDim() == r.GetDim() );

		for( GW_U32 i=0; i<r.GetDim(); ++i )
			r.SetData( i, a.GetData(i)+b.GetData(i) );
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  -
	*
	*  \param  v right hand statement
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        GW_VectorND  operator-(const GW_VectorND & v) const
	{
		GW_VectorND Res( this->GetDim() );
                 Minus( *this, v, Res );
		return Res;
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  Minus
	*
	*  \param  a left hand statement
	*  \param  b right hand statement
	*  \param  r result
	*  \author Gabriel Peyré 2001-09-29
	*
	*	Substract the two vectOMLs.
	*/ 
	/*------------------------------------------------------------------------------*/
        static void Minus(const GW_VectorND& a, const GW_VectorND& b, GW_VectorND& r)
	{
		GW_ASSERT( a.GetDim() == r.GetDim() );
		GW_ASSERT( b.GetDim() == r.GetDim() );

		for( GW_U32 i=0; i<r.GetDim(); ++i )
			r.SetData( i, a.GetData(i)-b.GetData(i) );
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name : operator  -
	*
	*  \return The opposite of the vector.
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        GW_VectorND  operator-() const
	{
		GW_VectorND Res( this->GetDim() );
                 UMinus( *this, Res );
		return Res;
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  UMinus
	*
	*  \param  a right hand statement
	*  \param  r result
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        static void UMinus(const GW_VectorND& a, GW_VectorND& r)
	{
		GW_ASSERT( a.GetDim() == r.GetDim() );

		for( GW_U32 i=0; i<r.GetDim(); ++i )
			r.SetData( i, -a.GetData(i) );
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  +=
	*
	*  \param  v rigth hand statement.
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        void  operator+=(const GW_VectorND & v)
	{	
                 Add( *this, v, *this );
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  -=
	*
	*  \param  v rigth hand statement.
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        void  operator-=(const GW_VectorND & v)
	{
                 Minus( *this, v, *this );
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  *
	*
	*  \param  v rigth hand statement.
	*	\return Dot product of the two vectOMLs.
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        GW_Float operator*(const GW_VectorND& v)
	{
		GW_ASSERT( this->GetDim() == v.GetDim() );

		GW_Float rDot=0;
		for( GW_U32 i=0; i<v.GetDim(); ++i )
			rDot += this->GetData(i)*v.GetData(i);

		return rDot;
	}


	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  NormInf
	*
	*  \return Infinite norm of the vector.
	*  \author Gabriel Peyré 2001-09-30
	*
	*	The infinte norm is the maximum of the absolute value of the coordonates.
	*/ 
	/*------------------------------------------------------------------------------*/
        GW_Float NormInf()
	{
		GW_Float rNorm = 0;
		GW_Float rVal = 0;
		for( GW_U32 i=0; i<this->GetDim(); ++i )
		{
			rVal = this->GetData(i);
			if( rVal>rNorm )
				rNorm = rVal;
			if( -rVal>rNorm )
				rNorm = -rVal;
		}

		return rNorm;
	}


	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  Norm2
	*
	*  \return Norm 2 of the vector.
	*  \author Gabriel Peyré 2001-09-30
	*
	*	The eclidian norm of the vector.
	*/ 
	/*------------------------------------------------------------------------------*/
        GW_Float Norm2()
	{
		GW_Float rNorm = 0;
		GW_Float rVal = 0;
		for( GW_U32 i=0; i<this->GetDim(); ++i )
		{
			rVal = this->GetData(i);
			rNorm += rVal*rVal;
		}

		return (GW_Float) sqrt( rNorm )/this->GetDim();
	}


	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  Norm1
	*
	*  \return Norm 1 of the vector.
	*  \author Gabriel Peyré 2001-09-30
	*
	*	The norm 1 is the sum of the absolute value of the coords.
	*/ 
	/*------------------------------------------------------------------------------*/
        GW_Float Norm1()
	{
		GW_Float rNorm = 0;
		for( GW_U32 i=0; i<this->GetDim(); ++i )
			rNorm += GW_ABS( this->GetData(i) );

		return rNorm/this->GetDim();
	}

	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  GetData
	*
	*  \param  i offset in the vector.
	*  \return the value of the given coord.
	*  \authOML Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        GW_Float GetData(GW_U32 i) const
	{
		GW_ASSERT( i<this->GetDim() );
		return (*this)[i];
	}


	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  SetData
	*
	*  \param  i offset in the vector.
	*  \param  rVal value to set to this coord.
	*  \authOML Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        void SetData(GW_U32 i, GW_Float rVal)
	{
		GW_ASSERT( i<this->GetDim() );
		(*this)[i] = rVal;
	}


	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  Randomize
	*
	*  \param  rMin=0 minimum value.
	*  \param  rMax=1 maximum value.
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        void  Randomize(GW_Float rMin = 0, GW_Float rMax = 1)
	{
		for( GW_U32 i=0; i<this->GetDim(); ++i )
			this->SetData( i, rMin + GW_RAND*(rMax-rMin) );
	}

	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  SetZero
	*
	*  \author Gabriel Peyré 2001-09-29
	*
	*	set the vector to zero.
	*/
	/*------------------------------------------------------------------------------*/
        void  SetZero()
	{
		this->SetValue(0.0);
	}
	/*------------------------------------------------------------------------------*/
	/** 
        * Name :  SetValue
	*
	*  \param  rVal the value to set to the vector.
	*  \author Gabriel Peyré 2001-09-29
	*/ 
	/*------------------------------------------------------------------------------*/
        void  SetValue(GW_Float rVal)
	{
		if( this->GetDim()>0 )
		{
			for( GW_U32 i=0; i<this->GetDim(); ++i )
				this->SetData(i,rVal);
		}
	}

	static void TestClass(std::ostream &s = cout)
	{
		TestClassHeader("GW_VectorND", s);

		TestClassFooter("GW_VectorND", s);
	}

private:

};


/*------------------------------------------------------------------------------*/
/** 
* Name :  operator
*
*  \param  f Rigth hand statement
*  \return the vector scaled.
*  \author Gabriel Peyré 2001-09-29
*/ 
/*------------------------------------------------------------------------------*/
inline 
GW_VectorND operator*(const GW_Float& f, GW_VectorND& v)
{
	GW_VectorND Res( v.GetDim() );
         GW_VectorND::Multiply( f, v, Res );
	return Res;
}

inline
std::ostream& operator<<(std::ostream &s, GW_VectorND& v)
{
	s << "GW_VectorND, size=" << v.GetDim() << " : " << endl;
	cout << "|";
	for( GW_U32 j=0; j<v.GetDim(); ++j )
	{
		s << v.GetData(j);
		if( j!=v.GetDim()-1 )
			s << " ";
	}
	s << "|";
	return s;
}



} // End namespace GW


#endif // _GW_VECTORND_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
