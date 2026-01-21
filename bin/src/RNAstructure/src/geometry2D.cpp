#include <math.h>
#include "geometry2D.h"

namespace geometry2D {
vec2D::vec2D() : x(0), y(0) { }
vec2D::vec2D(const vec2D& copy) : x(copy.x), y(copy.y) { }
vec2D::vec2D(const double x, const double y) : x(x), y(y) { }
vec2D::vec2D(const vec2D& from, const vec2D& to) : x(to.x-from.x), y(to.y-from.y) { }
vec2D vec2D::fromAngle(const double theta) {
	return vec2D(cos(theta), sin(theta));
}
vec2D vec2D::midpoint(const vec2D& a, const vec2D& b) {
	return vec2D ((a.x+b.x)/2, (a.y+b.y)/2);
}

vec2D& vec2D::set(const double x, const double y) { this->x=x;this->y=y; return *this;}
vec2D& vec2D::set(const vec2D& other) { x=other.x;y=other.y; return *this;}
vec2D& vec2D::add(const double dx, const double dy) { x+=dx;y+=dy; return *this;}
vec2D& vec2D::add(const vec2D& other) { x+=other.x;y+=other.y; return *this;}
vec2D& vec2D::subtr(const double dx, const double dy){ x-=dx;y-=dy; return *this;}
vec2D& vec2D::subtr(const vec2D& other){ x-=other.x;y-=other.y; return *this;}
vec2D& vec2D::scale(const double factor){ x*=factor;y*=factor; return *this;}
vec2D& vec2D::scale(const double scaleX, const double scaleY){ x*=scaleX;y*=scaleY; return *this;}
vec2D& vec2D::negate(){ x=-x;y=-y; return *this;}
//! Sets this vec2D's length to 1. 
//! If the vector's current length is 0, it will NOT be normalized. Instead it will be set to identically zero.
// Use normalize(double, double) to better handle the case when the vector's length is 0.
vec2D& vec2D::normalize() { return normalize(0, 0); }
vec2D& vec2D::normalize(const double x_if_Zero, const double y_if_Zero) {
    double mag = length();
    if (mag < 0.000000000001) {
        x = x_if_Zero;
        y = y_if_Zero;
    } else {
        x /= mag;
        y /= mag;
    }
    return *this;
}
	   
//! Sets the length (magnitude) of this vec2D.
//! If newLength is negative, the direction of the vec2D is also reversed.
//! This is equivalent to normalize() followed by scale(double)
//! 
//! Behavior may be unexpected if the current length is near zero.
vec2D& vec2D::setLength(double newLength) { scale(newLength / length()); return *this;}

//! Rotates this vec2D 90 degrees counter-clockwise.  
//! (x is replaced with -y and y is replaced with x).
//! The rotated vec2D is perpendicular (orthogonal) to its previous state.
//! This is equivalent to rotate(PI/2) or  rotate(-3*PI/2)
vec2D& vec2D::rotate90() { return set(-y, x); }

//! Rotates this vec2D 270 degrees counter-clockwise (or 90 degrees clockwise).
//! (x is replaced with y and y is replaced with -x).
//! The rotated vec2D is perpendicular (orthogonal) to its previous state.
//! This is equivalent to rotate(3*PI/2) or  rotate(-PI/2)
vec2D& vec2D::rotate270() { return set(y, -x); }
vec2D& vec2D::rotate(const double theta) {
    double cosT = cos(theta);
    double sinT = sin(theta);
    // x' = x cos theta - y sin theta
    // y' = x sin theta + y cos theta
    double xp = x;
    x = x * cosT - y * sinT;
    y = xp * sinT + y * cosT;
    return *this;
}

//! The magnitude (aka Euclidean length) of the vec2D. I.e: sqrt(x*x + y*y)
double vec2D::length() const { return sqrt(x*x + y*y); }
//! Returns the length-squared of this vector. I.e:   x*x + y*y
double vec2D::lengthSq() const { return x*x + y*y; }
//! Returns the Euclidean dot product of this vector with other. I.e:  A dot B = A.x * B.x + A.y * B.y
double vec2D::dot(const vec2D& other) const { return x * other.x + y * other.y; }
//! Returns the Euclidean dot product of this vector with (endPoint-referencePoint).
double vec2D::dot(const vec2D& endPoint, const vec2D& referencePoint) const { return x * (endPoint.x - referencePoint.x) + y * (endPoint.y - referencePoint.y); }
double vec2D::angle() const { return atan2(y, x); }

//! Returns a copy of this vec2D
vec2D vec2D::clone() const { return vec2D(x,y); }
//! Returns a new vec2D that is the midpoint between this and other.
vec2D vec2D::midpointTo(const vec2D& other) const { return midpoint(*this, other); }

//! Returns the location of result to be the point that would be reached by starting
//! at startPoint and traveling along the direction of this vec2D the distance specified by lengthAlongvec2D.
//! (The current length of this vec2D is irrelevant. I.e. It acts as if it were a unit vec2D.)
//! This is equivalent to {@code this.clone().setLength(lengthAlongvec2D).add(startPoint) }
vec2D vec2D::pointAlong(const vec2D& startPoint, const double lengthAlongvec2D) const { return clone().setLength(lengthAlongvec2D).add(startPoint); }
    //! Sets the location of result to be the point that would be reached by starting
    //!  at startPoint and traveling along the direction of this vec2D the specified multiple of this vec2D's length.
    //!  This is equivalent to {@code this.clone().scale(lengthMultiple).add(startPoint) }
    //!  For example  @{code new Vec2D(p1, p2).setPointAlongScaled(p1, 0.5, p3) } would set the
    //!  point p3 to be the midpoint between p1 and p2.
vec2D vec2D::pointAlongScaled(const vec2D& startPoint, double lengthMultiple) const {
    return clone().scale(lengthMultiple).add(startPoint);
}
//! Returns the projection of this vector along the specified direction.
//! This is the same as converting direction to a unit vector, U_dir, 
//! then multiplying U_dir by the dot product between this vector and U_dir.
vec2D vec2D::projectedOn(const vec2D& direction) const {
    return direction.clone().scale(dot(direction) / direction.lengthSq());
}

double vec2D::distanceTo(const vec2D& other) const { return distance(x, other.x, y, other.y); }
//! Gets the full angle from this vec2D to v2, sweeping counter-clockwise. The result is from -PI to PI. The order of the vec2Ds is important: {@code v1.angleTo(v2) == -v2.angleTo(v1)}  */
double vec2D::angleTo(const vec2D& other) const {
    // The dot product returns the smallest angle between the two vec2Ds,
    // so the result is always between 0 and PI.
    // But if we want to determine the full angle moving counter clockwise from v1 to v2,
    // we have to compare the angle to each vec2D.

	// Note: dot product = cos(a) * |v1||v2|
	// determinant = sin(a) * |v1||v2|  (same as length of cross product)
    double det = x*other.y-y*other.x; 
    return atan2(det, dot(other)); // i.e. atan(r sin(a), r cos(a)) where r=|v1||v2| is the magnitude an doesn't affect the angle.
}
	    
//! This calculates the angle from one vec2D to another, similar to angleTo(vec2D) except that the
//! output is intended for applications in which the direction of travel is indicated by the sign of the result
//! and the extent of travel is indicated by its absolute value. Thus instead of returning a value from -PI to PI, this
//! function returns a value from -2PI to 2PI.
//!
//! For example, want to draw an arc from 0 degrees to 270 degrees counter-clockwise. (I.e. the arc should span 3 quadrants.)
//! So v1 = (1, 0) and v2 = (0, -1)
//! But if {@link #angleTo(Vec2D)} were used, it would return -90 degrees, which would result in the arc being drawn CLOCKWISE through only the bottom-right quadrant.
//! This function would return 270 degrees (reverse = false) or -90 (reverse = true);
//!
//! As another example, assume v1 points to (-1, 0) and v2 points to (1, 0). This function would return 180 (reverse = false) or -180 (reverse = true);
//! As another example, assume v1 points to (1, 1) and v2 points to (-1, 1). This function would return 90 (reverse = false) or -270 (reverse = true);
//!
//! @param v2 the other vec2D
//! @param reverse whether the arc should be drawn counter-clockwise (reverse=false) or clockwise (reverse=true) from v1 to v2.
//!                         Alternatively reverse=true could be interpreted as "draw the arc counter-clockwise from v2 to v1 instead of from v1 to v2"
//!                If reverse is true, the output will be from -360 > angle >= 0
//!                If reverse is false, the output will be from 0 <= angle <= 360
double vec2D::arcAngleTo(const vec2D& other, bool reverse) const {
    double angle = angleTo(other);
    if (reverse && angle > 0)
        angle -= pi2; // return the equivalent "clockwise" (negative) angle
    else if (!reverse && angle < 0)
        angle += pi2; // return the equivalent "counter-clockwise" (positive) angle.
    return angle;
}


//-----------------------free functions-------------------------------------
//! Calculates the distance between two vec2Ds.
double distance(const vec2D& a, const vec2D& b) { return distance(a.x, a.y, b.x, b.y); }
//! Calculates the distance between two points with coordinates (x1, y1) and (x2, y2)
double distance(const double x1, const double y1, const double x2, const double y2){
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
//! Calculates the distance-squared between two points with coordinates (x1, y1) and (x2, y2)
double distanceSq(const double x1, const double y1, const double x2, const double y2) {
	return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}
//! Gets the smallest angle between two vec2Ds. 
//! The result is always between 0 and PI (inclusive). 
//! Based on the dot product A dot B = |A| * |B| * cos(theta)
double angleBetween(const vec2D& v1, const vec2D& v2) {
    return acos(v1.dot(v2) / sqrt(v1.lengthSq() * v2.lengthSq()));
}

//! Given any angle in radians, this returns an equivalent angle theta in the range [-pi, pi]
//! In general odd multiples of pi map to:  pi if they are positive and -pi if they are negative.
//! (Even multiples of pi, e.g. 2pi, map to 0.)
double wrapAngleToPi(double angle){
    return angle < 0 ?
        fmod((angle+pi), pi2) - pi : // returns -PI if angle is an exact multiple of PI
        fmod((angle-pi), pi2) + pi;  // returns PI if angle is an exact multiple of PI
}
//! Given any angle in radians, this returns an equivalent angle theta in the range [0, 2*pi] 
//! In general, multiples of 2*pi map to:  2*pi if they are positive and 0 if they are negative.
double wrapAngleTo2Pi(double angle){
    //return (angle = fmod(angle, pi2)) < 0 ? angle+pi2 : angle;
	return wrapAngleToPi(angle)+pi;
}
//! Given any angle in radians, this returns an equivalent angle theta in the range [-pi, pi)
//! Note that any multiple of pi maps to -pi because pi is an exclusive upper limit.
//! In order to have pi as an inclusive upper limit, use the wrapAngleToPi function.
double angleModPi(double angle){
        return fmod(angle, pi2);
}
//! Given any angle in radians, this returns an equivalent angle theta in the range [0, 2*pi) 
//! Note that multiples of 2*pi map to 0 because 2*pi is an exclusive upper limit.
//! In order to have 2*pi as an inclusive upper limit, use the wrapAngleTo2Pi function.
double angleMod2Pi(double angle){
    return fmod(angle, pi2)+pi;
}

vec2D& vec2D::operator=(const vec2D& other) { x=other.x; y=other.y; return *this; }
vec2D& vec2D::operator+=(const vec2D& other) { x+=other.x; y+=other.y; return *this; }
vec2D& vec2D::operator-=(const vec2D& other) { x-=other.x; y-=other.y; return *this; }
vec2D& vec2D::operator*=(const double scale) { x*=scale; y*=scale; return *this; }
vec2D& vec2D::operator/=(const double scale) { x/=scale; y/=scale; return *this; }

bool vec2D::operator==(const vec2D& other) { return x==other.x && y==other.y; }
bool vec2D::operator!=(const vec2D& other) { return x!=other.x || y!=other.y; }

//! Unary - operator: return the opposite (negative) vector
vec2D operator-(const vec2D& v) { return vec2D(-v.x, -v.y); }
vec2D operator+(const vec2D& v1, const vec2D& v2) { return vec2D(v1.x+v2.x, v1.y+v2.y); }
vec2D operator-(const vec2D& v1, const vec2D& v2) { return vec2D(v1.x-v2.x, v1.y-v2.y); }
vec2D operator*(const vec2D& v, const double scale) { return vec2D(v.x*scale, v.y*scale); }
vec2D operator*(const double scale, const vec2D& v) { return vec2D(v.x*scale, v.y*scale); }
vec2D operator/(const vec2D& v, const double scale) { return vec2D(v.x/scale, v.y/scale); }


//! Calculates the radius of a circle that has the given chord length and partial circumference.
//! Purpose: We have two points that lay on a circle of unknown radius. 
//! We know the straight-line distance between the points (aka the chord length) as well as the partial 
//! circumference from one point around the circle to the other point. But we do not know the remaining 
//! distance around the circle (on the other side of the points).
//! So the question is a simple one: what radius will yield a cicle with the given chord length and partial circumference?
//! This seems like an easy problem, but it turns out to be one that cannot be solved analytically.
//! Some definitions:
//!		S: length of the chord		 (known)    S = 'chard' = r * arcAngle
//!		c: circumference the circle  (unknown)  c = circWithoutArc + arcLength
//!		r: radius of the circle      (unknown)  r = c/(2*pi)  
//! This is used to calculate the optimal radius of an RNA loop that is cut by a bond on one (or more) side(s).
double calcRadiusFromChord(double circWithoutArc, double chord, double &arcAngle) {
	// Note: 1.110720735  is 2*PI()/4 * 1/(2*SIN(PI()/4)) 
	//       which is the factor to multiply chord by, assuming that the angle cut out by the chord is about 1/4th of the circle. 
	//       This is the largest chord:arcLength ratio that is reasonable assuming that a loop will contain at 
	//       least 3 nucleotides in addition to the two base-pair nucs.
	//       And bondLength should be roughly 2x or 3x the loop spacing.
	double circ = circWithoutArc + chord * 1.111;  // get an estimate of the hairpin loop circumference, assuming that chord cuts off about 1/4th of the circle.
	double r = circ/pi2;
	// circ was just an estimate because we did't yet know the arc-length across the bond (just the chord length). 
	arcAngle = 2 * asin(chord / 2 / r); // calculate angle from chord length  (using sin(T)=opp/hyp where opp=chord/2 and hyp=r)
	// here we recalculate to get a better estimate:
	circ = circWithoutArc+arcAngle*r; // get a BETTER estimate of the hairpin loop circumference. 
	r = circ/pi2;
	return r;
}

} // namespace geometry2D