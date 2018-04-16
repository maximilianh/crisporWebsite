/*
  Copyright by Matthias Hoechsmann (C) 2002-2004
  =====================================                                   
  You may use, copy and distribute this file freely as long as you
  - do not change the file,
  - leave this copyright notice in the file,
  - do not make any profit with the distribution of this file
  - give credit where credit is due
  You are not allowed to copy or distribute this file otherwise
  The commercial usage and distribution of this file is prohibited
  Please report bugs and suggestions to <mhoechsm@TechFak.Uni-Bielefeld.DE>
*/

#ifndef _ALGEBRA_H_
#define _ALGEBRA_H_

/** Algebra is a virtual class (interface) for algebra used by the Alignment template classes.
 *  This is in the spirit of the Algebraic Dynamic Programming (ADP) approach where a Dynamic Programming
 *  algorithm is seperated into a grammar and an algebra.
 */
template<class R, class L>
class Algebra
{
public:
    virtual R empty() const =0;								/**< Result for the empty tree alignment */
    virtual R replace(L a,R down, L b, R over) const =0;	/**< Result for the tree edit function 'replace' */    
    virtual R del(L a,R down, R over) const =0;				/**< Result for the tree edit function 'delete' */    
    virtual R insert(R down,L b,R over) const =0;			/**< Result for the tree edit function 'insert' */ 
    virtual R choice(R a,R b) const =0;						/**< The choice function. Commonly used functions are 'min' and 'max' */ 
    virtual R worst_score() const =0;						/**< The worst_score with respect to choice is specified by this function */

    virtual ~Algebra() {};
}; 

/** Extended Algebra for aligning RNA secondary structure trees where basepair replacements
 *  are considered as a single edit operation.
 */
template<class R, class L>
class RNA_Algebra : public Algebra<R,L>
{
public:
	/** Result for the replacement of a basepair */
	virtual R replacepair(L la, L lb, R down, L ra, L rb, R over) const =0;
};

/** Algebra for the tree edit model */
template<class R, class L>
class SZAlgebra
{
public:
    virtual R empty() const =0;								/**< Result for the empty tree alignment */
    virtual R replace(L a,R down, L b) const =0;			/**< Result for the tree edit function 'replace' */    
    virtual R del(L a,R down) const =0;						/**< Result for the tree edit function 'delete' */    
    virtual R insert(R down,L b) const =0;					/**< Result for the tree edit function 'insert' */ 
    virtual R choice(R a,R b) const =0;						/**< The choice function. Commonly used functions are 'min' and 'max' */ 
    virtual R worst_score() const =0;						/**< The worst_score with respect to choice is specified by this function */

    virtual ~SZAlgebra() {};
}; 

#endif





