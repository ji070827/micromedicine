"""Compatibility shim: provide torch_scatter API via PyG 2.8 builtins.

DiffDock 官方依赖 torch_scatter (C++ 扩展，Windows/Py3.13 无法编译)。
PyG 2.8 已内置 scatter 功能，此 shim 做等价映射。
"""
from torch_geometric.utils import scatter as _scatter


def scatter(src, index, dim=0, dim_size=None, reduce='sum'):
    return _scatter(src, index, dim=dim, dim_size=dim_size, reduce=reduce)


def scatter_mean(src, index, dim=0, dim_size=None):
    return _scatter(src, index, dim=dim, dim_size=dim_size, reduce='mean')


def scatter_sum(src, index, dim=0, dim_size=None):
    return _scatter(src, index, dim=dim, dim_size=dim_size, reduce='sum')


__all__ = ['scatter', 'scatter_mean', 'scatter_sum']