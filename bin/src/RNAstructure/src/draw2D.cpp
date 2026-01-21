#include "defines.h"
#include <math.h>
#include "structure.h"
#include "substructure.cpp"
#include "draw2D.h"

#define maxbranch 100
#define maxnopair 300  // maximum number of bases that can be forced single

namespace draw2D {

coordinates::coordinates(int size) {
  short int i;
  bases = size;
  x = new int[size + 1];
  y = new int[size + 1];
  num = new int*[size / 10 + 2];
  for (i = 0; i <= size / 10 + 1; i++) num[i] = new int[2];
}

coordinates::~coordinates() {
  short int i;
  delete[] x;
  delete[] y;
  for (i = 0; i <= bases / 10 + 1; i++) {
    delete[] num[i];
  }
  delete[] num;
}

void coordinates::set(int index, int x, int y) {
  this->x[index] = x;
  this->y[index] = y;
}
void coordinates::set(int index, double x, double y) {
  this->x[index] = (int)x;
  this->y[index] = (int)y;
}

// type is a function that determines the type of loop being entered by a 5'
//	nucleotide that is paired
//	type = 0 => pseudo-knot detected
//	type = 1 => hairpin
//	type = 2 => internal loop
//	type = 3 => multibranched loop
//	type = 4 => external loop
//	pseudoknots are excluded!!!
int type(int i, structure* ct, int structnum, int* intervene) {
  int j;
  // intervene counts the number of intervening helixes as the loops is followed
  //	back to its origin
  /*   if ((i==1)&&(ct->basepr[structnum][1]!=0)) {
 //this is an exterior loop with the first base paired
      j = ct->basepr[structnum][1];
    (*intervene)=1;
    while(true) {
      j++;
       if (j>ct->GetSequenceLength()) return 4;
       else if (ct->basepr[structnum][j]!=0) {
              (*intervene)++;
          j = ct->basepr[structnum][j];
       }
              }
 } */
  // else {
  j = i;
  (*intervene) = 0;
  while (true) {
    j++;
    // if (j >= ct->GetSequenceLength()) return 4;
    if (ct->GetPair(j, structnum) == i) {
      if ((*intervene) > 1) return 3;
      if ((*intervene) == 1)
        return 2;
      else
        return 1;
    }
    if (ct->GetPair(j, structnum) != 0) {
      if (++(*intervene) > ct->GetSequenceLength()) {
        cout << "Encountered Pseudoknot in loop type detection." << endl;
        return 0;
      }
      j = ct->GetPair(j, structnum);
    }
  }
  //   }
}
int type4(int i, structure* ct, int structnum, int* intervene) {
  int j;
  // intervene counts the number of intervening helixes as the loops is followed
  //	back to its origin
  if ((i == 1) && (ct->GetPair(1, structnum) != 0)) {
    // this is an exterior loop with the first base paired
    j = ct->GetPair(1, structnum);
    (*intervene) = 1;
    while (true) {
      j++;
      if (j > ct->GetSequenceLength())
        return 4;
      else if (ct->GetPair(j, structnum) != 0) {
        if (++(*intervene) > ct->GetSequenceLength()) {
          cout << "Encountered Pseudoknot in loop type detection." << endl;
          return 0;
        }
        j = ct->GetPair(j, structnum);
      }
    }
  } else {
    j = i;
    (*intervene) = 0;
    while (true) {
      j++;
      if (j >= ct->GetSequenceLength()) return 4;
      if (ct->GetPair(j, structnum) == i) {
        if ((*intervene) > 1) return 3;
        if ((*intervene) == 1)
          return 2;
        else
          return 1;
      }
      if (ct->GetPair(j, structnum) != 0) {
        if (++(*intervene) > ct->GetSequenceLength()) {
          cout << "Encountered Pseudoknot in loop type detection." << endl;
          return 0;
        }
        j = ct->GetPair(j, structnum);
      }
    }
  }
}

inline double distance(int x1, int y1, int x2, int y2) {
  return (sqrt((double)((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))));
}

inline double add(double a, double b) {  // function adds two angles
  if ((a + b) < 2 * pi)
    return (a + b);
  else
    return (a + b - 2 * pi);
}

// number calculates the coordinates that a number needs to be placed
void number(int i, coordinates* coord, double norm, int diameter) {
  coord->num[i / 10][1] = coord->y[i] + (int)(4 * diameter * cos(norm));
  // if (norm<pi) {
  coord->num[i / 10][0] = coord->x[i] + (int)(4 * diameter * sin(norm));
  //}
  // else {
  // if (coord->num[i/10][0]<coord->x[i]) {
  //	if (i<100)
  //  	coord->num[i/10][0] = coord->x[i] - 2*width;
  //	else if (i<1000)
  //  	coord->num[i/10][0] = coord->x[i] - 3*width;
  //	else
  // 	coord->num[i/10][0] = coord->x[i] - 4*width;
  //}
  return;
}

// place is used to output a graphical structure.
//	it puts the coordinates of the nucleotide positions
//	into the coord struct out.
//	These coordinates are output-type independent
void place(int structnum, structure* ct, coordinates* coord, int height,
           int width) {
  int diameter = (int)ceil(sqrt((double)((width)*width + (height)*height)));
  place(structnum, ct, coord, diameter, diameter * 0.4, diameter * 1 / 8.0f,
        diameter * 1 / 8.0f);
}
void place(int structnum, structure* ct, coordinates* coord, int diameter,
           float bondLength, float helixSpacing, float loopSpacing) {
  substructure sub;
  int size[maxbranch], looptype, x;
  int intervene, position[maxbranch];
  double ang, norm, angout[maxbranch], congestion[maxbranch], partial, highest,
      loopradius, xcenter, ycenter;
  double length1, length2, angle;
  double deltay, deltax;
  int i;  // the position of the nucleotide we are dealing with in the outer
  // loop.

  cout << "diameter:" << diameter << endl;
  cout << "bondLength:" << bondLength << endl;
  cout << "helixSpacing:" << helixSpacing << endl;
  cout << "loopSpacing:" << loopSpacing << endl;

  // ang = angle at which a helix leaves a loop
  // diameter = diameter of a single nucleotide character
  // norm = angle of a normal to the direction of helix propagation
  // size = size of a loop (number of nucs)
  // loopradius = radius of a loop
  // xcenter, ycenter = coordinates of a loop center

  // place a dummy nucleotide on the stack, this will put a space in the
  // external
  //	loop
  sub.putsub(0);
  coord->x[0] = 0;
  coord->y[0] = 0;
  // ct->basepr[structnum][0] = 0;

  while (sub.getsub(&i)) {
    if (ct->GetPair(i, structnum) != 0) {  // place a helix
      // find the angle at which the helix should proceed
      ang = pi * (i + ct->GetPair(i, structnum) - 1) / ct->GetSequenceLength();
      // angle of helix propagation
      norm = ang + pi / 2;

      if ((i % 10) == 0)  // place a number by nucleotide i
        number(i, coord, add(pi, norm), diameter);

      // place a's partner
      coord->x[ct->GetPair(i, structnum)] =
          coord->x[i] + (int)((diameter + bondLength) * sin(norm));
      coord->y[ct->GetPair(i, structnum)] =
          coord->y[i] + (int)((diameter + bondLength) * cos(norm));
      // now place the helix

      if (((ct->GetPair(i, structnum)) % 10) == 0)  // place a number
        number(ct->GetPair(i, structnum), coord, norm, diameter);
      i++;

      while (ct->GetPair(i, structnum) == (ct->GetPair(i - 1, structnum) - 1)) {
        coord->x[i] =
            coord->x[i - 1] + (int)((diameter + helixSpacing) * sin(ang));
        coord->y[i] =
            coord->y[i - 1] + (int)((diameter + helixSpacing) * cos(ang));

        if ((i % 10) == 0)  // place a number
          number(i, coord, add(norm, pi), diameter);

        // now place the base paired to i

        coord->x[ct->GetPair(i, structnum)] =
            coord->x[i] + (int)((diameter + bondLength) * sin(norm));
        coord->y[ct->GetPair(i, structnum)] =
            coord->y[i] + (int)((diameter + bondLength) * cos(norm));
        if ((ct->GetPair(i, structnum) % 10) == 0)  // place a number
          number(ct->GetPair(i, structnum), coord, norm, diameter);

        i++;
      }
    }
    // now we're in a loop, place the loop elements:
    //	distinguish between internal, hairpin, and multibranch loops

    ////////////////
    //   cout << "i - 1 ="<<i-1<<"\n";

    if (i == 0)
      looptype = type4(1, ct, structnum, &intervene) == 0 ? 0 : 4;
    else
      looptype = type(i - 1, ct, structnum, &intervene);
    //	type 0 => pseudo-knot detected
    //	type 1 => hairpin
    //	type 2 => internal loop
    //	type 3 => multibranched loop
    //	type 4 => external loop
    if (looptype == 0) {
      cout << "Exiting nucleobase placement due to pseudo-knot." << endl;
      return;
    }

    ///////////////////
    //   cout << "type = "<<looptype;
    //   cout << "\ninterevene = "<<intervene<<flush;
    // cin >> j;

    if (looptype == 1) {
      // hairpin, count the number of bases to determine the loop size
      size[0] = ct->GetPair(i - 1, structnum) - i;
      double angBond;
      loopradius = geometry2D::calcRadiusFromChord((size[0] + 1) * (diameter + loopSpacing),
                              diameter + bondLength, angBond);
      xcenter = coord->x[i - 1] + (diameter + bondLength) / 2 * sin(norm) +
                loopradius * sin(ang);
      ycenter = coord->y[i - 1] + (diameter + bondLength) / 2 * cos(norm) +
                loopradius * cos(ang);
      double angNuc =
          ang + pi + angBond / 2;  // start angle points at bond nucleotide
      double angStep =
          (diameter + loopSpacing) / loopradius;  // this should really be:
      // 2*pi*(diameter+loopSpacing)/circ  but
      // circ=loopradius*2*pi, so we can cancel out
      // the 2*pi and just leave loopradius.
      for (int j = 0; j < size[0]; j++) {  // place each nuc in turn
        angNuc += angStep;  // rotate partially around the circle.
        coord->x[i + j] = (int)(xcenter + loopradius * sin(angNuc));
        coord->y[i + j] = (int)(ycenter + loopradius * cos(angNuc));
        if ((i + j) % 10 == 0)  // place a number
          number(i + j, coord, angNuc, diameter);
      }
    } else if (looptype == 2) {
      // internal loop, count the number of bases on each side
      // where is the next helix?
      x = i;
      while (ct->GetPair(x, structnum) == 0) x++;
      size[0] = x - i;
      size[1] = ct->GetPair(i - 1, structnum) - ct->GetPair(x, structnum) - 1;
      // cout << "size0 = "<<size[0]<<"\n" << "size1 = "<<size[1]<<"\n";

      // calculate which side is more congested
      // angout = angle of outgoing helix
      if (ang < pi)
        angout[0] = ang + pi;
      else
        angout[0] = ang - pi;

      angout[1] =
          pi * (x + ct->GetPair(x, structnum) - 1) / ct->GetSequenceLength();
      // cout << "ang = "<<ang/pi<<"\n";
      // cout << "angout[0] = "<<angout[0]/pi<<"\n";
      // cout << "angout[1] = "<<angout[1]/pi<<"\n";

      congestion[0] = size[0] + 1 / (fabs(angout[1] - angout[0]));
      congestion[1] = size[1] + 1 / (fabs(angout[0] - angout[1]));
      // cout << "c0 = "<<congestion[0]<<"\n";
      // cout << "c1 = "<<congestion[1]<<"\n";

      // double angDiff =
      // the congested side will determine the seperation of helixes
      if (congestion[0] > congestion[1])
        loopradius = (size[0] + 2) * (diameter + loopSpacing) /
                     pi2;  //(5*fabs(angout[1] - angout[0]))
      else
        loopradius = (size[1] + 2) * (diameter + loopSpacing) /
                     pi2;  ///(5*fabs(angout[0] - angout[1]))

      // place the elements:
      xcenter = coord->x[i - 1] + (diameter + bondLength) / 2 * sin(norm) +
                loopradius * sin(ang);
      ycenter = coord->y[i - 1] + (diameter + bondLength) / 2 * cos(norm) +
                loopradius * cos(ang);

      if (angout[0] < angout[1])
        partial = fabs(angout[1] - angout[0]);
      else
        partial = fabs(2 * pi + angout[1] - angout[0]);

      for (int j = 1; j <= size[0]; j++) {
        double angNuc = angout[0] + partial * j / (1 + size[0]);
        coord->x[i + j - 1] = (int)(xcenter + loopradius * sin(angNuc));
        coord->y[i + j - 1] = (int)(ycenter + loopradius * cos(angNuc));

        if ((i + j - 1) % 10 == 0)  // place a number
          number(i + j - 1, coord, angNuc, diameter);
      }
      for (int j = 1; j <= size[1]; j++) {
        if (angout[1] < angout[0])
          partial = fabs(angout[0] - angout[1]);
        else
          partial = fabs(2 * pi - angout[1] + angout[0]);

        ////////
        // cout << "j = "<<j;
        // cout << "  partial angle = "<< (angout[1]+
        // fabs(angout[0]-angout[1])*j/(1+size[1]))/pi<<"\n";

        coord->x[ct->GetPair(x, structnum) + j] =
            (int)(xcenter +
                  loopradius *
                      sin(angout[1] +
                          partial * ((double)j) / (1 + (double)size[1])));
        coord->y[ct->GetPair(x, structnum) + j] =
            (int)(ycenter +
                  loopradius *
                      cos(angout[1] +
                          partial * ((double)j) / (1 + (double)size[1])));

        if ((ct->GetPair(x, structnum) + j) % 10 == 0)  // place a number
          number(ct->GetPair(x, structnum) + j, coord,
                 (angout[1] + partial * (j) / (1 + size[1])), diameter);
      }
      // Place the outgoing helix on the stack
      sub.putsub(x);

      ////////
      // cout << "putting x = "<<x;
      // cin >> j;

      coord->x[x] = (int)(xcenter - (diameter + bondLength) / 2 * sin(norm) +
                          loopradius * sin(angout[1]));
      coord->y[x] = (int)(ycenter - (diameter + bondLength) / 2 * cos(norm) +
                          loopradius * cos(angout[1]));
    } else {  // external or multibranch loop
      if (looptype == 4)
        position[0] = 0;  // case of external loop
      else
        position[0] = ct->GetPair(i - 1, structnum);
      for (int j = 1; j <= intervene; j++) {
        position[j] = ct->GetPair(position[j - 1], structnum);
        (position[j])++;
        while (ct->GetPair(position[j], structnum) == 0) (position[j])++;
        size[j - 1] = position[j] - ct->GetPair(position[j - 1], structnum) - 1;
      }

      if (looptype == 4)
        size[intervene] = ct->GetSequenceLength() -
                          ct->GetPair(position[intervene], structnum);
      else
        size[intervene] =
            position[0] - ct->GetPair(position[intervene], structnum) - 1;

      // calculate which side is more congested
      // angout = angle of outgoing helix
      if (looptype == 4) {
        angout[0] = 0;
        ang = pi;
        norm = 3 * pi / 2;
      } else {
        if (ang < pi)
          angout[0] = ang + pi;
        else
          angout[0] = ang - pi;
      }

      angout[intervene + 1] = angout[0];
      for (int j = 1; j <= intervene; j++) {
        angout[j] = pi *
                    (position[j] + ct->GetPair(position[j], structnum) - 1) /
                    ct->GetSequenceLength();
        congestion[j - 1] =
            (size[j - 1] + 1) / (fabs(angout[j] - angout[j - 1]));
      }
      congestion[intervene] =
          (size[intervene] + 1) / (fabs(angout[0] - angout[intervene]));
      highest = congestion[0];
      x = 0;
      for (int j = 1; j <= intervene; j++) {
        if (congestion[j] > highest) {
          highest = congestion[j];
          x = j;  // most congested side
        }
      }

      // the congested side will determine the seperation of helixes
      loopradius = (size[x] + 2) * diameter * (2 * pi) /
                   (5 * fabs(angout[x + 1] - angout[x]));

      if (looptype == 4) norm = 3 * pi / 2;

      if (looptype == 4) {
        xcenter =
            loopradius * sin(ang) + (diameter + bondLength) / 2 * sin(norm);
        ycenter =
            loopradius * cos(ang) + (diameter + bondLength) / 2 * cos(norm);
      } else {
        xcenter = coord->x[i - 1] + loopradius * sin(ang) +
                  (diameter + bondLength) / 2 * sin(norm);
        ycenter = coord->y[i - 1] + loopradius * cos(ang) +
                  (diameter + bondLength) / 2 * cos(norm);
      }

      for (x = 0; x <= intervene; x++) {
        for (int j = 1; j <= size[x]; j++) {
          if (angout[x] < angout[x + 1])
            partial = fabs(angout[x + 1] - angout[x]);
          else
            partial = fabs(2 * pi + angout[x + 1] - angout[x]);

          coord->x[ct->GetPair(position[x], structnum) + j] =
              (int)(xcenter +
                    loopradius * sin(angout[x] +
                                     partial * ((double)(j + 1)) /
                                         ((double)(2 + size[x]))));

          ////////
          // cout << "partial angle = "<< (angout[0]+
          // fabs(angout[1]-angout[0])*j/(1+size[0]))<<"\n";

          coord->y[ct->GetPair(position[x], structnum) + j] =
              (int)(ycenter +
                    loopradius * cos(angout[x] +
                                     partial * ((double)(j + 1)) /
                                         ((double)(2 + size[x]))));

          if ((ct->GetPair(position[x], structnum) + j) % 10 ==
              0)  // place a number
            number(ct->GetPair(position[x], structnum) + j, coord,
                   (angout[x] +
                    partial * ((double)(j + 1)) / ((double)(2 + size[x]))),
                   diameter);
        }
      }

      // position the helixes:
      for (int j = 1; j <= intervene; j++) {
        coord->x[position[j]] =
            (int)(xcenter +
                  loopradius * sin(angout[j]));  //- 0.8*diameter*sin(norm)
        coord->y[position[j]] =
            (int)(ycenter +
                  loopradius * cos(angout[j]));  //- 0.8*diameter*cos(norm)
      }

      // Do a first order correction on the placement of nucs in a multibranch
      //	loop
      // only do correction if the loops are not congested at all:

      if (looptype == 3 &&
          ((congestion[0] * diameter < (loopradius / 2)) &&
           (congestion[intervene] * diameter < (loopradius / 2)))) {
        // find the distance to from the entering helix to the first an last
        //	exiting helixes

        // distance from entering helix to first exiting helix
        length1 = distance(coord->x[ct->GetPair(position[0], structnum)],
                           coord->y[ct->GetPair(position[0], structnum)],
                           coord->x[position[1]], coord->y[position[1]]);

        // distance from entering helix to last exiting helix
        // first place the the 3' nuc in the last exiting helix
        angle = pi * (position[intervene] +
                      ct->GetPair(position[intervene], structnum) - 1) /
                ct->GetSequenceLength();
        // angle of helix propagation
        norm = angle + pi / 2;

        // place the 3' nuc
        coord->x[ct->GetPair(position[intervene], structnum)] =
            coord->x[position[intervene - 1]] +
            (int)((diameter + bondLength) * sin(norm));
        coord->y[ct->GetPair(position[intervene], structnum)] =
            coord->y[position[intervene - 1]] +
            (int)((diameter + bondLength) * cos(norm));

        length2 =
            distance(coord->x[position[0]], coord->y[position[0]],
                     coord->x[ct->GetPair(position[intervene], structnum)],
                     coord->y[ct->GetPair(position[intervene], structnum)]);

        // now decide which side is most congested:
        congestion[1] = size[0] / length1;

        congestion[2] = size[intervene] / length2;

        if (congestion[1] > congestion[2]) {  // first side more congested
          // find angle that nucs need to be laid out on
          if (size[0] > 0) {
            ////debugging from abs to fabs caused overlap
            angle = fabs(angout[1] - ang) * loopradius /
                    (diameter * (double)size[0]);
            for (int j = 1; j <= size[0]; j++) {
              coord->x[ct->GetPair(position[0], structnum) + j] =
                  coord->x[ct->GetPair(position[0], structnum) + j - 1] +
                  (int)(diameter * sin(ang - angle));
              coord->y[ct->GetPair(position[0], structnum) + j] =
                  coord->y[ct->GetPair(position[0], structnum) + j - 1] +
                  (int)(diameter * cos(ang - angle));
              if ((ct->GetPair(position[0], structnum) + j) % 10 == 0) {
                if (angle + ang > pi / 2) {
                  norm = angle + ang - pi / 2;
                } else
                  norm = angle + ang + 3 * pi / 2;
                number(ct->GetPair(position[0], structnum) + j, coord, norm,
                       diameter);
              }
            }
          }
          // how much closer can things be moved?
          length1 =
              distance(coord->x[ct->GetPair(position[0], structnum) + size[0]],
                       coord->y[ct->GetPair(position[0], structnum) + size[0]],
                       coord->x[position[1]], coord->y[position[1]]) -
              diameter;

          if (size[intervene] > 0) {
            angle = fabs(ang - angout[intervene]) * loopradius /
                    (diameter * size[intervene]);
            for (int j = size[intervene]; j > 0; j--) {
              coord->x[ct->GetPair(position[intervene], structnum) + j] =
                  coord
                      ->x[ct->GetPair(position[intervene], structnum) + j + 1] +
                  (int)(diameter * sin(angle + ang));
              coord->y[ct->GetPair(position[intervene], structnum) + j] =
                  coord
                      ->y[ct->GetPair(position[intervene], structnum) + j + 1] +
                  (int)(diameter * cos(angle + ang));

              if ((ct->GetPair(position[intervene], structnum) + j) % 10 == 0) {
                norm = add(pi / 2, angle + ang);
                number(ct->GetPair(position[intervene], structnum) + j, coord,
                       norm, diameter);
              }
            }
          }
          // else {
          //	length1=distance(coord->x[position[0]+size[1]],
          //	coord->y[position[0]+size[1]],coord->x[position[1]],
          //	coord->y[position[1]]) - diameter;
          //}
        } else {  // back side more congested:
          // find angle that nucs need to be laid out on
          if (size[intervene] > 0) {
            angle = fabs(ang - angout[intervene]) * loopradius /
                    (diameter * size[intervene]);
            for (int j = size[intervene]; j > 0; j--) {
              coord->x[ct->GetPair(position[intervene], structnum) + j] =
                  coord
                      ->x[ct->GetPair(position[intervene], structnum) + j + 1] +
                  (int)(diameter * sin(angle + ang));
              coord->y[ct->GetPair(position[intervene], structnum) + j] =
                  coord
                      ->y[ct->GetPair(position[intervene], structnum) + j + 1] +
                  (int)(diameter * cos(angle + ang));
              if ((ct->GetPair(position[intervene], structnum) + j) % 10 == 0) {
                norm = add(pi / 2, angle + ang);
                number(ct->GetPair(position[intervene], structnum) + j, coord,
                       norm, diameter);
              }
            }
          }
          // how much closer can things be moved?
          length1 =
              distance(
                  coord->x[ct->GetPair(position[intervene], structnum) + 1],
                  coord->y[ct->GetPair(position[intervene], structnum) + 1],
                  coord->x[ct->GetPair(position[intervene], structnum)],
                  coord->y[ct->GetPair(position[intervene], structnum)]) -
              diameter;

          if (size[0] > 0) {
            angle = fabs(angout[1] - ang) * loopradius / (diameter * size[0]);
            for (int j = 1; j <= size[0]; j++) {
              coord->x[ct->GetPair(position[0], structnum) + j] =
                  coord->x[ct->GetPair(position[0], structnum) + j - 1] +
                  (int)(diameter * sin(ang - angle));
              coord->y[ct->GetPair(position[0], structnum) + j] =
                  coord->y[ct->GetPair(position[0], structnum) + j - 1] +
                  (int)(diameter * cos(ang - angle));
              if ((ct->GetPair(position[0], structnum) + j) % 10 == 0) {
                if (angle + ang > pi / 2) {
                  norm = angle + ang - pi / 2;
                } else
                  norm = angle + ang + 3 * pi / 2;
                number(ct->GetPair(position[0], structnum) + j, coord, norm,
                       diameter);
              }
            }
          }
          // else {
          // length1 = distance(coord->x[position[intervene-1]+1],
          //	coord->y[position[intervene-1]+1],coord->x[position[intervene-1]],
          //	coord->y[intervene-1]) - diameter;
          //}
        }

        // move all the elements closer
        for (x = 1; x <= intervene; x++) {
          coord->x[position[x]] =
              coord->x[position[x]] - int(length1 * (sin(ang)));
          coord->y[position[x]] =
              coord->y[position[x]] - int(length1 * (cos(ang)));
        }
        for (x = 1; x <= intervene - 1; x++) {
          for (int j = 1; j <= size[x]; j++) {
            coord->x[ct->GetPair(position[x], structnum) + j] =
                coord->x[ct->GetPair(position[x], structnum) + j] -
                (int)(length1 * sin(ang));
            coord->y[ct->GetPair(position[x], structnum) + j] =
                coord->y[ct->GetPair(position[x], structnum) + j] -
                (int)(length1 * cos(ang));

            if ((ct->GetPair(position[x], structnum) + j) % 10 ==
                0) {  // move a number
              coord->num[(ct->GetPair(position[x], structnum) + j) / 10][0] =
                  coord
                      ->num[(ct->GetPair(position[x], structnum) + j) / 10][0] -
                  int(length1 * (sin(ang)));
              coord->num[(ct->GetPair(position[x], structnum) + j) / 10][1] =
                  coord
                      ->num[(ct->GetPair(position[x], structnum) + j) / 10][1] -
                  int(length1 * (cos(ang)));
            }
          }
        }

        // place the nucs between the external-most helix and the two adjacent
        // helices in a line connecting their position
        // 3' to external helix, # nucs = size[0]
        // 5' to external helix, # nucs = size [intervene]
        // Added 8/2/04 by D.H.M.

        // start 3'
        deltax = coord->x[position[1]] -
                 coord->x[ct->GetPair(position[0], structnum)];
        deltay = coord->y[position[1]] -
                 coord->y[ct->GetPair(position[0], structnum)];
        deltax = deltax / ((double)size[0] + 1.0);
        deltay = deltay / ((double)size[0] + 1.0);
        for (int j = 1; j <= size[0]; j++) {
          coord->x[ct->GetPair(position[0], structnum) + j] =
              coord->x[ct->GetPair(position[0], structnum)] +
              (int)(((double)j) * deltax);
          coord->y[ct->GetPair(position[0], structnum) + j] =
              coord->y[ct->GetPair(position[0], structnum)] +
              (int)(((double)j) * deltay);
        }

        // now do the 5' nucs
        /*deltax=coord->x[position[intervene]]-coord->x[ct->basepr[structnum][position[intervene-1]]];
        deltay=coord->y[position[intervene]]-coord->y[ct->basepr[structnum][position[intervene-1]]];
        deltax=deltax/((double) size[0]+1.0);
        deltay=deltay/((double) size[0]+1.0);
        for (int j=1;j<=size[intervene];j++) {
                coord->x[ct->basepr[structnum][position[intervene-1]]+j]=coord->x[ct->basepr[structnum][position[intervene-1]]]+(int)
        (((double) j)*deltax );
                coord->y[ct->basepr[structnum][position[intervene-1]]+j]=coord->y[ct->basepr[structnum][position[intervene-1]]]+(int)
        (((double) j)*deltay );
        }*/
      }

      // Place the outgoing helixes on the stack
      for (int j = 1; j <= intervene; j++) {
        // cout << "here?\n";
        sub.putsub(position[j]);
      }
    }
  }

  // As a final step, remove labels that will superimpose on nucleotides
  for (i = 10; i <= ct->GetSequenceLength(); i = i + 10) {
    for (int j = 1; j <= ct->GetSequenceLength(); j++) {
      if (coord->x[j] + diameter > coord->num[i / 10][0] &&
          coord->x[j] < coord->num[i / 10][0] + diameter) {
        if (coord->y[j] + diameter > coord->num[i / 10][1] &&
            coord->y[j] < coord->num[i / 10][1] + diameter) {
          // conflict in location:
          coord->num[i / 10][0] = 0;
          coord->num[i / 10][1] = 0;
        }
      }
    }
  }
}

// place the coordinates for a structure with pseudoknots.  This uses a circle.
void placepk(structure* ct, coordinates* coord, int height, int width) {
  double diagonal;
  double radius;
  int i;

  // This is needed as legacy code for support of place(), above.
  coord->x[0] = 0;
  coord->y[0] = 0;

  // First determine the size of a diagonal for a character:
  diagonal = sqrt(((double)height) * ((double)width));
  radius = 0.2 * diagonal * ((double)ct->GetSequenceLength());

  // place the nucs:
  for (i = 1; i <= ct->GetSequenceLength(); ++i) {
    coord->x[i] =
        -radius * sin(2 * pi * ((double)i) / ((double)ct->GetSequenceLength()));
    coord->y[i] =
        -radius * cos(2 * pi * ((double)i) / ((double)ct->GetSequenceLength()));
  }

  // Place the numbers:
  for (i = 10; i <= ct->GetSequenceLength(); i += 10) {
    coord->num[i / 10][0] =
        -(radius + 5 * diagonal) *
        sin(2 * pi * ((double)i) / ((double)ct->GetSequenceLength()));
    coord->num[i / 10][1] =
        -(radius + 5 * diagonal) *
        cos(2 * pi * ((double)i) / ((double)ct->GetSequenceLength()));
  }
}

void sortxy(coordinates* coord, bool counter, int height, int width) {
  // Translate structure to origin.
  // Draw Clockwise if required.
  int i, x, y;
  int diameter = (int)sqrt((double)(width)*width + (height)*height);
  if (!counter) {  // translate coordinates to make a clockwise structure
    // go through every coordinate in out
    for (i = 1; i <= coord->bases; i++) {
      coord->x[i] = -coord->x[i];
    }
    for (i = 10; i <= coord->bases; i = i + 10) {
      coord->num[i / 10][0] = -coord->num[i / 10][0];
    }
  }

  // sort the coordinates for the lowest x and y values
  x = coord->x[0];
  y = coord->y[0];

  for (i = 1; i <= coord->bases; i++) {
    if (coord->x[i] < x) x = coord->x[i];
    if ((i % 10) == 0 &&
        !(coord->num[i / 10][0] == 0 && coord->num[i / 10][1] == 0)) {
      // Note, 0,0 is flag for number coords that conflict with BP coords
      if (coord->num[i / 10][0] < x) x = coord->num[i / 10][0];
      if (coord->num[i / 10][1] < y) y = coord->num[i / 10][1];
    }
    if (coord->y[i] < y) y = coord->y[i];
  }

  // now place the coordinates into out, translating so that they are 1 diameter
  //	from 0,0
  x = x - diameter;
  y = y - diameter;

  for (i = 1; i <= coord->bases; i++) {
    coord->x[i] = coord->x[i] - x;
    coord->y[i] = coord->y[i] - y;
  }

  for (i = 10; i <= coord->bases; i = i + 10) {
    if (!(coord->num[i / 10][0] == 0 && coord->num[i / 10][1] == 0)) {
      coord->num[i / 10][0] = coord->num[i / 10][0] - x;
    }
  }
}

} // namespace draw2D