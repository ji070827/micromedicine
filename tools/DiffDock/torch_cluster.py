"""Compatibility shim: provide torch_cluster API via PyG 2.8 builtins.

DiffDock 官方依赖 torch_cluster (C++ 扩展) 的 radius/radius_graph/knn_graph。
PyG 2.8 内置这些功能，此 shim 做等价映射。
"""
from torch_geometric.nn import radius_graph, knn_graph
from torch_geometric.nn.pool import radius


__all__ = ['radius', 'radius_graph', 'knn_graph']