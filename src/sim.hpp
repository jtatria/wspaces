/*
 * Copyright (C) 2017 José Tomás Atria <jtatria at gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SIMM_
#define SIMM_ 1

#include "matrix.hpp"

inline double simm_add( const Vec vi, const Vec vj, const int i, const int j ) {
    double ii  = vi.coeff( i ), ij = vj.coeff( i ), ji = vi.coeff( j ), jj = vj.coeff( j );
    double d   = ( ( ii > 0 && ji > 0 ? ii + ji : 0 ) + ( ij > 0 && jj > 0 ? ij + jj : 0 ) );
    double num = (
        ( vi.array().sign() * vj.array().sign() ).array() * ( vi + vj ).array()
    ).sum() - d;
    double den = ( vi.sum() - ( ( ii > 0 ? ii : 0 ) + ( ii > 0 ? ii : 0 ) ) );
    return num / den;
}

inline double simm_dw( const Vec vi, const Vec vj, const int i, const int j ) {
    double ii  = vi.coeff( i ), ij = vj.coeff( i ), ji = vi.coeff( j ), jj = vj.coeff( j );
    double d   = ( std::min( ii, ij ) + std::min( ji, jj ) );
    double num = vi.cwiseMin( vj ).sum() - d;
    double den = vi.sum() - ( ii + ji );
    return num / den;
}

inline double simm_cos( const Vec vi, const Vec vj, const int i, const int j ) {
    double num = ( vi.dot( vj ) );
    double den =  ( vi.norm() * vj.norm() );
    return num / den;
}

#endif // SIMM_
