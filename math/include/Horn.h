#ifndef _HORN_H_
#define _HORN_H_

#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>
#include <fstream>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>

#include "math/include/dataobject.h"

#ifndef PI
	#ifndef M_PI
		#define PI 3.1415926535897932384626433832795
	#else
		#define PI M_PI
	#endif
#endif

/** 	\class	Horn
 *	\brief	Stellt die Mittel zur Berechnung der Transformation zwischen zwei System
 *
 * 	Die Klasse Horn dient der Suche nach der besten
 * 	Transformation zwischen zwei KoordinatenSystemen (Left und Right)
 * 	Das linke Koordinatensystem in unserem Fall ist der KUKA, das rechte System das Tracking-System
 * 	der Kalibrierungsvorgang erfolgt nach dem Paper: "Closed-form solution of absolute orientation using unit quaternions"
 */

namespace Math {

	class Horn {
	private:
		void Alloc(void);
		/**
		*	\brief	Berechnet die Schwerpunkte der beiden Punktwolken
		*
		*	Die Methode berechnet die Schwerpunkte der Punktwolken, die 
		*	in _dataPointsLeft und _dataPointsRight stehen.
		*
		*	@param 	rightCentroid	Schwerpunkt vom rechten System wird hier gespeichert
		*	@param	leftCentroid	Schwerpunkt vom linken System wird hier gespeichert
		*/
		void findCentroids(gsl_vector* rightCentroid, gsl_vector* leftCentroid);

		/**
		*	\brief	relativiert die Punkte bez&uuml;glich ihres Schwerpunktes
		*
		*	Von den Punkten in den Matrizen _dataPointsLeft und _dataPointsRight
		*	werden die Schwerpunkte abgezogen
		*
		*/
		void relativizeDataPoints(gsl_vector* rightCentroid, gsl_vector* leftCentroid);

		/**
		*	\brief	Berechnung einer 3x3-Matrix M, die die Informationen aller Punkte enth&auml;t
		*
		*	Die Matrix enth&auml;lt die Informationen aller Punkte und wird nach folgender Vorschrift berechnet:
		*	\f[	M = \sum_{i=1}^{\_n} l_i^T r_i \f]
		*	wobei \f$ l_i \f$ der i-te Punkt der linken und \f$ r_i \f$ der i-te Punkt der rechten relativierten Punktwolken sei.
		*
		*	@param	matrixM	Speicher f&uuml;r die Matrix M
		*/
		void construct3x3Matrix(gsl_matrix* matrixM);

		/**
		*	\brief Berechnung der Matrix N aus dem Paper
		*
		*	Die Matrix N hat folgenden Aufbau. <br>
		*	Sei \f$ M = (m)_{ij} \f$ mit \f$ i=1,2,3 \f$
		*		\f[ N = \left ( \begin{array}{cccc} m_{11}+m_{22}+m_{33} & m_{23}-m_{32}        & m_{31}-m_{13}         & m_{12}-m{21} \\ 
		*					  	m_{23}-m_{32}         & m_{11}-m_{22}-m_{33} & m_{12}+m_{21}         & m_{31}+m_{13} \\ 
		*						m_{31}-m_{13}         & m_{12}+m_{21}        & -m_{11}+m_{22}-m_{33} & m_{23}+m_{32} \\ 
		*						m_{12}-m{21}          & m_{31}+m_{13}        & m_{23}+m_{32}         & -m_{11}-m_{22}+m_{33} \\ \end{array} \right ) 
		*	\f]
		*/
		void construct4x4Matrix(gsl_matrix* matrixN, gsl_matrix* matrixM);

		/**
		*	\brief 	Berechnet den Eigenvektor zum gr&ouml;&szlig;ten Eigenwert von N
		*
		*	Der Eigenwert zum zum gr&ouml;&szlig;ten Eigenwert von N ist das Quaternion,
		*	dass den Rotationsanteil der Transformation vom linken zum rechten System beschreibt.
		*
		*	@param	matrixM	Die 3x3Matrix M 
		*	@param	Speicher f&uuml;r die Matrix N
		*/
		void findEigenvector(gsl_matrix* matrixN);

		/**
		*	\brief	Berechnung der Rotationsmatrix aus dem Eigenvektor
		*
		*	Die Rotationsmatrix wird aus dem Eigenvektor (Einheits-Quaternion berechnet).
		*	Dadruch ist sie "quasi" orthonormal.
		*/
		void buildRotationMatrix();

		/**
		*	\brief 	Berechnung des Skalierungsfaktors
		*
		*	Es wird der symmetrische Skalierungsfaktor aus den Punktpaaren berechnet.
		*	Es muss nichts weiter bekannt sein. 
		*/
		void calculateScale();

		/**
		*	\brief	Berechnung der Translationsanteile
		*
		*	Die Methode berechnet die beiden Translationsanteile der Transformation (links->rechts und umgekehrt)
		*	Die Translation t wird mittels der Schwerpunkte berechnet, z.B. links->rechts 
		*	\f[ t =  C_r - R*C_l \f] 
		*	mit \f$ C_r = \f$ rightCentroid und \f$ C_l = \f$ leftCentroid
		*
		*	@param	leftCentroid	Schwerpunkt der linken Punktwolke
		*	@param	rightCentroid	Schwerpunkt der rechten Punktwolke
		*/
		void calculateTranslation(gsl_vector* rightCentroid, gsl_vector* leftCentroid);

		/**
		*	\brief	berechnet Maximum zweier Werte
		*
		*
		*
		*/
		double dmax(double a, double b);

		gsl_matrix*	_dataPointsLeft;  ///< Die Spalten der Matrix sind die Punkte des linken Systems (sp&auml;ter relativiert)
		gsl_matrix*	_dataPointsRight; ///< Die Spalten der Matrix sind die Punkte des rechten Systems (sp&auml;ter relativiert)
		gsl_matrix*	_dataPointsLeftOld; ///< Die Spalten der Matrix sind die Punkte des linken Systems
		gsl_matrix*	_dataPointsRightOld; ///< Die Spalten der Matrix sind die Punkte des rechten Systems
		gsl_matrix*	_matrixLeftToRight; ///< Die Rotationsmatrix für die Transformation vom linken zum rechten System
		double		_maxPosEigenvalue; ///<  Der gr&ouml;&szlig;te Eigenwert der berechneten 4x4-Matrix (siehe Paper, construct4x4Matrix)
		gsl_vector*	_eigenvector;	///< Der Eigenvektor zum gr&ouml;&szlig;ten Eigenwert, entspricht Quaternion der Transformation
		gsl_vector*	_translationLtoR; ///< Der Translationsanteil vom linken zum rechten System
		double		_scale;	///< Der symmetrische Skalierungfaktor der Transformation, Indikator f&uuml;r die G&uuml;te (1 sehr gut), da keine Skalierung vorliegen sollte
		int		_n; ///< Anzahl der Punktpaare

		gsl_vector* _leftCentroid;
		gsl_vector* _rightCentroid;
		gsl_matrix* _matrixM;
		gsl_matrix* _matrixN;

		gsl_vector* _eigenvalues;
		gsl_matrix* _eigenvectors;
		gsl_eigen_symmv_workspace* _work;

		gsl_vector* _rotatedLeft;
		gsl_vector* _rotatedRight;
		gsl_vector* _errorVector;

		bool _isCalculated;

		bool _needToFreePoints;

	public:
		/**
		 *	\brief Konstruktor der Klasse
		 *	
		 *	
		 *	@param	n	Anzahl der Punktpaare
		 *	@param	xMax	x-Radius der Punktwolke	
		 *	@param	yMax 	y-Radius der Punktwolke
		 *	@param	zMax	z-Radius der Punktwolke
		 *	@param	trackerServerIP	IP-Adresse des Servers des Trackingsystems
		 *	@param	trackerServerPort Port des Servers des Trackingsystems
		 *	@param	tracker	zu verwendender Marker
  		 *	
 		 */
		Horn(int n, gsl_matrix* dataPointsLeft, gsl_matrix* dataPointsRight);
		Horn(std::vector<tDataObject>& left, std::vector<tDataObject>& right);
		Horn(std::vector<Matrix>& left, std::vector<Matrix>& right);
		Horn(std::vector<DoubleVector>& left, std::vector<Matrix>& right);
		Horn(std::vector<Matrix>& left, std::vector<DoubleVector>& right);
		Horn(std::vector<DoubleVector>& left, std::vector<DoubleVector>& right);
		Horn(void);

		~Horn();

		void Init(int n, gsl_matrix* dataPointsLeft, gsl_matrix* dataPointsRight);

		/**
		*	\brief Brechnet die beste Transformation von zwei Koordinationsystemen
		*
		*	Die Berechnung der besten Transformation zwischen 2 Systemen erfolgt mittels
		*	einer Methode nach dem Paper "Closed-form solution of absolute orientation using unit quaternions".
		*	 
		*/
		void calculateTransformation(void);
//		double calculateTransformation(int type);

		/**
			 *	\brief	Berechnung des Fehlers 
  		 *	
  		 *	Berechnet den Fehler mittels der Summe der Quadrate.
		 *
		 */
		double errorSumOfSquares(std::vector<std::vector<double> > *errors = NULL);

		/**
			 *	\brief Der maximale Fehler 
  		 *
  		 *	Berechnet den maximalen Fehler, der mit der errechneten Transformation bei
		 *	den Punktpaare besteht.
		 *
		 */
		double maxError(void);

		bool getQuaternion(double quat[4]);
		bool getTranslation(double trans[3]);
		bool getScale(double &scale);
		bool getTransform(double quat[4], double trans[3], double &scale);
		bool getTransform(double matrixRowWise[12]);
		bool getTransform(Matrix &matrix);
		bool getMeanError(double &error);
		bool getMaxError(double &error);
	};

}

#endif
