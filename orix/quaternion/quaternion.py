# -*- coding: utf-8 -*-
# Copyright 2018-2022 the orix developers
#
# This file is part of orix.
#
# orix is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# orix is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with orix.  If not, see <http://www.gnu.org/licenses/>.

import warnings

import dask.array as da
from dask.diagnostics import ProgressBar
import numpy as np
import quaternion

from orix.base import check, Object3d
from orix.vector import Miller, Vector3d


def check_quaternion(obj):
    return check(obj, Quaternion)


class Quaternion(Object3d):
    r"""Basic quaternion object.

    Quaternions support the following mathematical operations:

    - Unary negation.
    - Inversion.
    - Multiplication with other quaternions and vectors.

    Quaternion-quaternion multiplication for two quaternions
    :math:`q_1 = (a_1, b_1, c_1, d_1)`
    and :math:`q_2 = (a_2, b_2, c_2, d_2)`
    with :math:`q_3 = (a_3, b_3, c_3, d_3) = q_1 * q_2` follows as:

    .. math::
       a_3 = (a_1 * a_2 - b_1 * b_2 - c_1 * c_2 - d_1 * d_2)

       b_3 = (a_1 * b_2 + b_1 * a_2 + c_1 * d_2 - d_1 * c_2)

       c_3 = (a_1 * c_2 - b_1 * d_2 + c_1 * a_2 + d_1 * b_2)

       d_3 = (a_1 * d_2 + b_1 * c_2 - c_1 * b_2 + d_1 * a_2)

    Quaternion-vector multiplication with a three-dimensional vector
    :math:`v = (x, y, z)` calculates a rotated vector
    :math:`v' = q * v * q^{-1}` and follows as:

    .. math::
       v'_x = x(a^2 + b^2 - c^2 - d^2) + 2(z(a * c + b * d) + y(b * c - a * d))

       v'_y = y(a^2 - b^2 + c^2 - d^2) + 2(x(a * d + b * c) + z(c * d - a * b))

       v'_z = z(a^2 - b^2 - c^2 + d^2) + 2(y(a * b + c * d) + x(b * d - a * c))

    Attributes
    ----------
    data : numpy.ndarray
        The numpy array containing the quaternion data.
    a, b, c, d : numpy.ndarray
        The individual elements of each vector.
    conj : Quaternion
        The conjugate of this quaternion :math:`q^* = a - bi - cj - dk`.
    """

    dim = 4

    @property
    def a(self):
        return self.data[..., 0]

    @a.setter
    def a(self, value):
        self.data[..., 0] = value

    @property
    def b(self):
        return self.data[..., 1]

    @b.setter
    def b(self, value):
        self.data[..., 1] = value

    @property
    def c(self):
        return self.data[..., 2]

    @c.setter
    def c(self, value):
        self.data[..., 2] = value

    @property
    def d(self):
        return self.data[..., 3]

    @d.setter
    def d(self, value):
        self.data[..., 3] = value

    @property
    def antipodal(self):
        return self.__class__(np.stack([self.data, -self.data], axis=0))

    @property
    def conj(self):
        q = quaternion.from_float_array(self.data).conj()
        return Quaternion(quaternion.as_float_array(q))

    def __invert__(self):
        return self.__class__(self.conj.data / (self.norm**2)[..., np.newaxis])

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            q1 = quaternion.from_float_array(self.data)
            q2 = quaternion.from_float_array(other.data)
            return other.__class__(quaternion.as_float_array(q1 * q2))
        elif isinstance(other, Vector3d):
            # check broadcast shape is correct before calculation, as
            # quaternion.rotat_vectors will perform outer product
            # this keeps current __mul__ broadcast behaviour
            q1 = quaternion.from_float_array(self.data)
            v = quaternion.as_vector_part(
                (q1 * quaternion.from_vector_part(other.data)) * ~q1
            )
            if isinstance(other, Miller):
                m = other.__class__(xyz=v, phase=other.phase)
                m.coordinate_format = other.coordinate_format
                return m
            else:
                return other.__class__(v)
        return NotImplemented

    def __neg__(self):
        return self.__class__(-self.data)

    @classmethod
    def triple_cross(cls, q1, q2, q3):
        """Pointwise cross product of three quaternions.

        Parameters
        ----------
        q1, q2, q3 : orix.quaternion.Quaternion
            Three quaternions for which to find the "triple cross".

        Returns
        -------
        q : orix.quaternion.Quaternion
        """
        q1a, q1b, q1c, q1d = q1.a, q1.b, q1.c, q1.d
        q2a, q2b, q2c, q2d = q2.a, q2.b, q2.c, q2.d
        q3a, q3b, q3c, q3d = q3.a, q3.b, q3.c, q3.d
        a = (
            +q1b * q2c * q3d
            - q1b * q3c * q2d
            - q2b * q1c * q3d
            + q2b * q3c * q1d
            + q3b * q1c * q2d
            - q3b * q2c * q1d
        )
        b = (
            +q1a * q3c * q2d
            - q1a * q2c * q3d
            + q2a * q1c * q3d
            - q2a * q3c * q1d
            - q3a * q1c * q2d
            + q3a * q2c * q1d
        )
        c = (
            +q1a * q2b * q3d
            - q1a * q3b * q2d
            - q2a * q1b * q3d
            + q2a * q3b * q1d
            + q3a * q1b * q2d
            - q3a * q2b * q1d
        )
        d = (
            +q1a * q3b * q2c
            - q1a * q2b * q3c
            + q2a * q1b * q3c
            - q2a * q3b * q1c
            - q3a * q1b * q2c
            + q3a * q2b * q1c
        )
        q = cls(np.vstack((a, b, c, d)).T)
        return q

    def dot(self, other):
        """Dot product of this quaternion and the other as a
        numpy.ndarray.
        """
        return np.sum(self.data * other.data, axis=-1)

    def dot_outer(self, other):
        """Outer dot product of this quaternion and the other as a
        numpy.ndarray.
        """
        dots = np.tensordot(self.data, other.data, axes=(-1, -1))
        return dots

    def mean(self):
        """Calculates the mean quaternion with unitary weights.

        Notes
        -----
        The method used here corresponds to Equation (13) in
        https://arc.aiaa.org/doi/pdf/10.2514/1.28949.
        """
        q = self.flatten().data.T
        qq = q.dot(q.T)
        w, v = np.linalg.eig(qq)
        w_max = np.argmax(w)
        return self.__class__(v[:, w_max])

    def outer(self, other, lazy=False, chunk_size=20, progressbar=True):
        """Compute the outer product of this quaternion and the other
        quaternion or vector.

        Parameters
        ----------
        other : orix.quaternion.Quaternion or orix.vector.Vector3d
        lazy : bool, optional
            Whether to computer this computation using Dask. This option
            can be used to reduce memory usage when working with large
            arrays. Default is False.
        chunk_size : int, optional
            When using `lazy` computation, `chunk_size` represents the
            number of objects per axis for each input to include in each
            iteration of the computation. Default is 20.
        progressbar : bool, optional
            Whether to show a progressbar during computation if `lazy`
            is True. Default is True.

        Returns
        -------
        orix.quaternion.Quaternion or orix.vector.Vector3d
        """

        if isinstance(other, Quaternion):
            if lazy:
                darr = self._outer_dask(other, chunk_size=chunk_size)
                arr = np.empty(darr.shape)
                if progressbar:
                    with ProgressBar():
                        da.store(darr, arr)
                else:
                    da.store(darr, arr)
            else:
                q1 = quaternion.from_float_array(self.data)
                q2 = quaternion.from_float_array(other.data)
                # np.outer works with flattened array
                q = np.outer(q1, q2).reshape(q1.shape + q2.shape)
                arr = quaternion.as_float_array(q)
            return other.__class__(arr)
        elif isinstance(other, Vector3d):
            if lazy:
                darr = self._outer_dask(other, chunk_size=chunk_size)
                arr = np.empty(darr.shape)
                if progressbar:
                    with ProgressBar():
                        da.store(darr, arr)
                else:
                    da.store(darr, arr)
            else:
                q = quaternion.from_float_array(self.data)
                arr = quaternion.rotate_vectors(q, other.data)
            if isinstance(other, Miller):
                m = other.__class__(xyz=arr, phase=other.phase)
                m.coordinate_format = other.coordinate_format
                return m
            else:
                return other.__class__(arr)
        else:
            raise NotImplementedError(
                "This operation is currently not avaliable in orix, please use outer "
                "with `other` of type `Quaternion` or `Vector3d`"
            )

    def _outer_dask(self, other, chunk_size=20):
        """Compute the product of every quaternion in this instance to
        every quaternion or vector in another instance, returned as a
        Dask array.

        For quaternion-quaternion multiplication, this is known as the
        Hamilton product.

        Parameters
        ----------
        other : orix.quaternion.Quaternion or orix.vector.Vector3d
        chunk_size : int, optional
            Number of objects per axis for each input to include in each
            iteration of the computation. Default is 20.

        Returns
        -------
        dask.array.Array

        Notes
        -----
        For quaternion-quaternion multiplication, to create a new
        quaternion from the returned array `arr`, do
        `q = Quaternion(arr.compute())`. Likewise for quaternion-vector
        multiplication, to create a new vector from the returned array
        do `v = Vector3d(arr.compute())`.
        """
        if not isinstance(other, (Quaternion, Vector3d)):
            raise TypeError("Other must be Quaternion or Vector3d.")

        ndim1 = len(self.shape)
        ndim2 = len(other.shape)

        # Set chunk sizes
        chunks1 = (chunk_size,) * ndim1 + (-1,)
        chunks2 = (chunk_size,) * ndim2 + (-1,)

        # Dask has no dask.multiply.outer(), use dask.array.einsum
        # Summation subscripts
        str1 = "abcdefghijklm"[:ndim1]  # Max. object dimension of 13
        str2 = "nopqrstuvwxyz"[:ndim2]
        sum_over = f"...{str1},{str2}...->{str1 + str2}"

        # Get quaternion parameters as dask arrays to be computed later
        q1 = da.from_array(self.data, chunks=chunks1)
        a1, b1, c1, d1 = q1[..., 0], q1[..., 1], q1[..., 2], q1[..., 3]

        # We silence dask's einsum performance warnings for "small"
        # chunk sizes, since using the chunk sizes suggested floods
        # memory
        warnings.filterwarnings("ignore", category=da.PerformanceWarning)

        if isinstance(other, Quaternion):
            q2 = da.from_array(other.data, chunks=chunks2)
            a2, b2, c2, d2 = q2[..., 0], q2[..., 1], q2[..., 2], q2[..., 3]
            # fmt: off
            a = (
                + da.einsum(sum_over, a1, a2)
                - da.einsum(sum_over, b1, b2)
                - da.einsum(sum_over, c1, c2)
                - da.einsum(sum_over, d1, d2)
            )
            b = (
                + da.einsum(sum_over, b1, a2)
                + da.einsum(sum_over, a1, b2)
                - da.einsum(sum_over, d1, c2)
                + da.einsum(sum_over, c1, d2)
            )
            c = (
                + da.einsum(sum_over, c1, a2)
                + da.einsum(sum_over, d1, b2)
                + da.einsum(sum_over, a1, c2)
                - da.einsum(sum_over, b1, d2)
            )
            d = (
                + da.einsum(sum_over, d1, a2)
                - da.einsum(sum_over, c1, b2)
                + da.einsum(sum_over, b1, c2)
                + da.einsum(sum_over, a1, d2)
            )
            # fmt: on
            out = da.stack((a, b, c, d), axis=-1)
        else:  # Vector3d
            v2 = da.from_array(other.data, chunks=chunks2)
            x2, y2, z2 = v2[..., 0], v2[..., 1], v2[..., 2]
            # fmt: off
            x = (
                + da.einsum(sum_over, a1**2 + b1**2 - c1**2 - d1**2, x2)
                + da.einsum(sum_over, a1 * c1 + b1 * d1, z2) * 2
                + da.einsum(sum_over, b1 * c1 - a1 * d1, y2) * 2
            )
            y = (
                + da.einsum(sum_over, a1**2 - b1**2 + c1**2 - d1**2, y2)
                + da.einsum(sum_over, a1 * d1 + b1 * c1, x2) * 2
                + da.einsum(sum_over, c1 * d1 - a1 * b1, z2) * 2
            )
            z = (
                + da.einsum(sum_over, a1**2 - b1**2 - c1**2 + d1**2, z2)
                + da.einsum(sum_over, a1 * b1 + c1 * d1, y2) * 2
                + da.einsum(sum_over, b1 * d1 - a1 * c1, x2) * 2
            )
            # fmt: on
            out = da.stack((x, y, z), axis=-1)

        new_chunks = tuple(chunks1[:-1]) + tuple(chunks2[:-1]) + (-1,)
        return out.rechunk(new_chunks)
