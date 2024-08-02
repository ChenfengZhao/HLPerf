import simpy
from utils import *
import timeit
import argparse

def read_nod_src(env, fifo_out_nod_src_list, nidx_begin, nidx_end, L, II, L_mem, deg_list):
    """read nod_src from memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_out_nod_src_list : list of simpy.store
        The list of output FIFO of neighbor index range (or degree) of each node . length is 6
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory read access
    deg_list : list of int

    """
    print("Simulation starts at %d" % env.now)

    for n in range(nidx_begin, nidx_end):
        # yield env.timeout(L_mem) # memory latency
        # yield env.timeout(15) # memory latency

        if n == nidx_begin:
            yield env.timeout(L+L_mem)
        else:
            yield env.timeout(II)

        for fifo_idx in range(len(fifo_out_nod_src_list)):
            yield fifo_out_nod_src_list[fifo_idx].put(deg_list[n])

    print("read_nod_src ends at %d" % env.now)

def read_edge_src(env, fifo_in_nod_src, fifo_out_tmp_src, nidx_begin, nidx_end, L, II, L_mem):
    """read edge_src from memory
    
    Instead of putting actual edge index into the FIFO, we put fake values

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_tmp_src : simpy.store
        Output FIFO of neibhbor indices of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory read access
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        # yield env.timeout(L_mem) # memory latency
        # yield env.timeout(15) # memory latency

        for e in range(deg):
            if e == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            yield fifo_out_tmp_src.put(1)

    print("read_edge_src ends at %d" % env.now)

def read_ft_in_tar(env, fifo_in_nod_src, fifo_out_ft_in_ah, fifo_out_ft_in_eh, nidx_begin, nidx_end, L, II, L_mem):
    """read ft_in for layer A and E from memory (target features)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_ft_in_ah : simpy.store
        Output FIFO of ah feature of each node
    fifo_out_ft_in_eh : simpy.store
        Output FIFO of feature of eh of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory read access
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        # yield env.timeout(L_mem) # memory latency

        for i in range(2):
            yield fifo_out_ft_in_ah.put(1)

            if i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            if deg != 0:
                yield fifo_out_ft_in_eh.put(1)
                # pass

    print("read_ft_in_tar ends at %d" % env.now)

def update_tar_ah(env, fifo_in_ft_in_ah, fifo_out_rst_ah_p1, nidx_begin, nidx_end, L_list, II_list, ft_size = 32, d = 8):
    """execute linear A layer (linear projection: in_feats -> out_feats)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_in_ah : simpy.store
        Input FIFO of ah feature of each node
    fifo_out_rst_ah_p1 : simpy.store
        Output FIFO of result feature (phase 1) of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    ft_size : int
        feature size of each node and edge
    d : int
        The segment size of Vector-matrix multiplication
    """

    for n in range(nidx_begin, nidx_end):

        for i in range(2):
            _ = yield fifo_in_ft_in_ah.get()
            if i == 0:
                yield env.timeout(L_list[0])
            else:
                yield env.timeout(II_list[0])

        yield env.timeout(L_list[1] + (ft_size-1) * II_list[1])

        for kd in range(d):
            if kd == 0:
                yield env.timeout(L_list[2])
            else:
                yield env.timeout(II_list[2])
            
            yield fifo_out_rst_ah_p1.put(1)
    
    print("update_tar_ah ends at %d" % env.now)

def update_tar_ah_sum(env, fifo_in_rst_ah_p1, fifo_out_rst_ah, nidx_begin, nidx_end, L_list, II_list, ft_size = 32, d = 8):
    """execute linear A layer - phase 2 (linear projection: in_feats -> out_feats)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_ah_p1 : simpy.store
        Input FIFO of result feature (phase 1) of each node
    fifo_out_rst_ah : simpy.store
        Output FIFO of result feature (phase 2) of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    ft_size : int, optional
        feature size of each node and edge, by default 32
    d : int, optional
        The segment size of Vector-matrix multiplication, by default 8
    """

    for n in range(nidx_begin, nidx_end):

        for kd in range(d):
            _ = yield fifo_in_rst_ah_p1.get()
            if kd == 0:
                yield env.timeout(L_list[0])
            else:
                yield env.timeout(II_list[0])
        
        # yield env.timeout(L_list[1] + (int(ft_size * d / 8) - 1) * II_list[1])
        yield env.timeout(L_list[1] + (ft_size * d / 8 - 1) * II_list[1])

        for i in range(2):
            if i == 0:
                yield env.timeout(L_list[2])
            else:
                yield env.timeout(II_list[2])
        
            yield fifo_out_rst_ah.put(1)

    print("update_tar_ah_sum ends at %d" % env.now)

def update_tar_eh(env, fifo_in_ft_in_eh, fifo_in_nod_src, fifo_out_rst_eh_p1, nidx_begin, nidx_end, L_list, II_list, ft_size = 32, d = 8):
    """execute linear E layer (linear projection: in_feats -> out_feats)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_in_eh : simpy.store
        Input FIFO of eh feature of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_eh_p1 : simpy.store
        Output FIFO of result eh feature (phase 1) of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    ft_size : int
        feature size of each node and edge
    d : int
        The segment size of Vector-matrix multiplication
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        if deg != 0:
            for i in range(2):
                _ = yield fifo_in_ft_in_eh.get()
                if i == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            yield env.timeout(L_list[1] + (ft_size-1) * II_list[1])

            for kd in range(d):
                if kd == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
                
                yield fifo_out_rst_eh_p1.put(1)

    print("update_tar_eh ends at %d" % env.now)

def update_tar_eh_sum(env, fifo_in_rst_eh_p1, fifo_in_nod_src, fifo_out_rst_eh, nidx_begin, nidx_end, L_list, II_list, ft_size = 32, d = 8):
    """execute linear E layer - phase 2 (linear projection: in_feats -> out_feats)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_eh_p1 : simpy.store
        Input FIFO of result feature (phase 1) of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_eh : _type_
        Output FIFO of result feature (phase 2) of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    ft_size : int, optional
        feature size of each node and edge, by default 32
    d : int, optional
        The segment size of Vector-matrix multiplication, by default 8
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        if deg != 0:
            for kd in range(d):
                _ = yield fifo_in_rst_eh_p1.get()
                if kd == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            yield env.timeout(L_list[1] + (int(ft_size * d / 8) - 1) * II_list[1])

            for i in range(2):
                if i == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
            
                yield fifo_out_rst_eh.put(1)

    print("update_tar_eh_sum ends at %d" % env.now)

def read_feat_in_agg(env, fifo_in_nod_src, fifo_in_tmp_src, fifo_out_ft_in_dh, fifo_out_ft_in_bh, nidx_begin, nidx_end, L, II, L_mem):
    """read aggregated features from memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_in_tmp_src : simpy.store
        Input FIFO of neighbor indices of each node
    fifo_out_ft_in_dh : simpy.store
        Output FIFO of input dh feature to be aggregated
    fifo_out_ft_in_bh : simpy.store
        Output FIFO of bh feature to be aggregated
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory read access
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for e in range(deg):
            _ = yield fifo_in_tmp_src.get()
            yield env.timeout(3)

            # yield env.timeout(L_mem) # memory latency

            for i in range(2):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)
                
                yield fifo_out_ft_in_dh.put(1)
                yield fifo_out_ft_in_bh.put(1)
    
    print("read_feat_in_agg ends at %d" % env.now)

def update_dh(env, fifo_in_ft_in_dh, fifo_in_nod_src, fifo_out_rst_dh_p1, nidx_begin, nidx_end, L_list, II_list, ft_size = 32, d = 8):
    """execute linear D layer (linear projection: in_feats -> out_feats)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_in_dh : simpy.store
        Input FIFO of dh feature
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_dh_p1 : simpy.store
        Output FIFO of result dh feature (phase 1)
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    ft_size : int
        feature size of each node and edge
    d : int
        The segment size of Vector-matrix multiplication
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for e in range(deg):

            for i in range(2):
                _ = yield fifo_in_ft_in_dh.get()
                if i == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            yield env.timeout(L_list[1] + (ft_size-1) * II_list[1])

            for kd in range(d):
                if kd == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
                
                yield fifo_out_rst_dh_p1.put(1)

    print("update_dh ends at %d" % env.now)

def update_dh_sum(env, fifo_in_rst_dh_p1, fifo_in_nod_src, fifo_out_rst_dh, nidx_begin, nidx_end, L_list, II_list, ft_size = 32, d = 8):
    """execute linear D layer - phase 2 (linear projection: in_feats -> out_feats)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_dh_p1 : simpy.store
        Input FIFO of result dh feature (phase 1) of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_dh : simpy.store
        Output FIFO of result dh feature (phase 2) of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    ft_size : int, optional
        feature size of each node and edge, by default 32
    d : int, optional
        The segment size of Vector-matrix multiplication, by default 8
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for e in range(deg):
            for kd in range(d):
                _ = yield fifo_in_rst_dh_p1.get()
                if kd == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            yield env.timeout(L_list[1] + (int(ft_size * d / 8) - 1) * II_list[1])

            for i in range(2):
                if i == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
            
                yield fifo_out_rst_dh.put(1)

    print("update_dh_sum ends at %d" % env.now)

def update_bh(env, fifo_in_ft_in_bh, fifo_in_nod_src, fifo_out_rst_bh_p1, nidx_begin, nidx_end, L_list, II_list, ft_size = 32, d = 8):
    """execute linear B layer (linear projection: in_feats -> out_feats)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_in_bh : simpy.store
        Input FIFO of bh feature
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_bh_p1 : simpy.store
        Output FIFO of result feature (phase 1)
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    ft_size : int
        feature size of each node and edge
    d : int
        The segment size of Vector-matrix multiplication
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for e in range(deg):
            for i in range(2):
                _ = yield fifo_in_ft_in_bh.get()
                if i == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            yield env.timeout(L_list[1] + (ft_size-1) * II_list[1])

            for kd in range(d):
                if kd == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
                
                yield fifo_out_rst_bh_p1.put(1)

    print("update_bh ends at %d" % env.now)

def update_bh_sum(env, fifo_in_rst_bh_p1, fifo_in_nod_src, fifo_out_rst_bh, nidx_begin, nidx_end, L_list, II_list, ft_size = 32, d = 8):
    """execute linear B layer - phase 2 (linear projection: in_feats -> out_feats)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_bh_p1 : simpy.store
        Input FIFO of result feature (phase 1) of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree)
    fifo_out_rst_bh : simpy.store
        Output FIFO of result feature (phase 2)
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    ft_size : int, optional
        feature size of each node and edge, by default 32
    d : int, optional
        The segment size of Vector-matrix multiplication, by default 8
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for e in range(deg):
            for kd in range(d):
                _ = yield fifo_in_rst_bh_p1.get()
                if kd == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            yield env.timeout(L_list[1] + (int(ft_size * d / 8) - 1) * II_list[1])

            for i in range(2):
                if i == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
            
                yield fifo_out_rst_bh.put(1)

    print("update_bh_sum ends at %d" % env.now)

def read_e_idx(env, fifo_in_nod_src, fifo_out_e_idx_list,nidx_begin, nidx_end, L, II, L_mem):
    """read e_idx from memeory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_e_idx_list : list of simpy.store
        The list of output FIFO of e_idx of each edge
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory read access
    """
    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        # yield env.timeout(L_mem) # memory latency

        for e in range(deg):
            if e == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            for i in range(len(fifo_out_e_idx_list)):
                yield fifo_out_e_idx_list[i].put(1)


def read_ft_in_edge_agg(env, fifo_in_nod_src, fifo_out_e_idx, fifo_out_ft_in_ce, nidx_begin, nidx_end, L, II, L_mem):
    """read aggregated edge features from memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_e_idx : simpy.store
        Output FIFO of e_idx of each edge
    fifo_out_ft_in_ce : simpy.store
        Output FIFO of feature of ce of each edge
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory read access
    """
    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for e in range(deg):
            _ = yield fifo_out_e_idx.get()
            # yield env.timeout(2)
            
            # yield env.timeout(L_mem) # memory latency
            for i in range(2):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)
                
                yield fifo_out_ft_in_ce.put(1)

    print("read_ft_in_edge_agg ends at %d" % env.now)

def update_ce(env, fifo_in_ft_in_ce, fifo_in_nod_src, fifo_out_rst_ce_p1, nidx_begin, nidx_end, L_list, II_list, ft_size = 32, d = 8):
    """execute linear C layer: edge (linear projection: in_feats -> out_feats)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_in_ce : simpy.store
        Input FIFO of ce feature of each edge
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_ce_p1 : simpy.store
        Output FIFO of result feature (phase 1) of each edge
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    ft_size : int
        feature size of each node and edge
    d : int
        The segment size of Vector-matrix multiplication
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for e in range(deg):
            for i in range(2):
                _ = yield fifo_in_ft_in_ce.get()
                if i == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])

            yield env.timeout(L_list[1] + (ft_size-1) * II_list[1])

            for kd in range(d):
                if kd == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
                
                yield fifo_out_rst_ce_p1.put(1)
    
    print("update_ce ends at %d" % env.now)

def update_ce_sum(env, fifo_in_rst_ce_p1, fifo_in_nod_src, fifo_out_rst_ce, nidx_begin, nidx_end, L_list, II_list, ft_size = 32, d = 8):
    """execute linear C layer - phase 2 (linear projection: in_feats -> out_feats)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_ce_p1 : simpy.store
        Input FIFO of result ce feature (phase 1) of each edge
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node
    fifo_out_rst_ce : simpy.store
        Output FIFO of result ce feature (phase 2) of each edge
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    ft_size : int, optional
        feature size of each node and edge, by default 32
    d : int, optional
        The segment size of Vector-matrix multiplication, by default 8
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for e in range(deg):
            yield env.timeout(3)
            for kd in range(d):
                _ = yield fifo_in_rst_ce_p1.get()
                if kd == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            yield env.timeout(L_list[1] + (int(ft_size * d / 8) - 1) * II_list[1])

            for i in range(2):
                if i == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
            
                yield fifo_out_rst_ce.put(1)
    
    print("update_ce_sum ends at %d" % env.now)

def comp_rst_e_sum_sigma(env, fifo_in_nod_src, fifo_in_rst_bh, fifo_in_rst_dh, fifo_in_rst_eh, fifo_in_rst_ce, fifo_out_rst_e, fifo_out_rst_sum_sigma, fifo_out_rst_sum_sigma_h, nidx_begin, nidx_end, L_list, II_list):
    """calculate rst_e, sum_sigma and sum_sigma_h

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node
    fifo_in_rst_bh : simpy.store
        Input FIFO of bh result feature
    fifo_in_rst_dh : simpy.store
        Input FIFO of dh result feature
    fifo_in_rst_eh : simpy.store
        Input FIFO of eh result feature
    fifo_in_rst_ce : simpy.store
        Input FIFO of ce result feature
    fifo_out_rst_e : simpy.store
        Output FIFO of result feature of each edge
    fifo_out_rst_sum_sigma : simpy.store
        Output FIFO of result feature of sum_sigma
    fifo_out_rst_sum_sigma_h : simpy.store
        Output FIFO of result feature of sum_sigma_h
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        if deg != 0:
            for i in range(2):
                _ = yield fifo_in_rst_eh.get()
                if i == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])

            for e in range(deg):
                for i in range(2):
                    if e == 0 and i == 0:
                        yield env.timeout(L_list[1])
                    else:
                        yield env.timeout(II_list[1])
                    
                    _ = yield fifo_in_rst_dh.get()
                    _ = yield fifo_in_rst_ce.get()
                    yield fifo_out_rst_e.put(1)
                    _ = yield fifo_in_rst_bh.get()
            
            for i in range(2):
                if i == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
                
                yield fifo_out_rst_sum_sigma.put(1)
                yield fifo_out_rst_sum_sigma_h.put(1)

    print("comp_rst_e_sum_sigma ends at %d" % env.now)

def write_rst_e_mem(env, fifo_in_e_idx, fifo_in_nod_src, fifo_in_rst_e, nidx_begin, nidx_end, L, II, L_mem):
    """write rst_e to memory for each edges

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_e_idx : simpy.store
        Input FIFO of 
    fifo_in_nod_src : simpy.store
        _description_
    fifo_in_rst_e : simpy.store
        _description_
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory write access
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for e in range(deg):
            _ = yield fifo_in_e_idx.get()
            yield env.timeout(1)

            # yield env.timeout(L_mem) # memory latency

            for i in range(2):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)
                
                _ = yield fifo_in_rst_e.get()

    print("write_rst_e_mem ends at %d" % env.now)

def comp_rst(env, fifo_in_rst_ah, fifo_in_nod_src, fifo_in_rst_sum_sigma_h, fifo_in_rst_sum_sigma, fifo_out_rst_h, nidx_begin, nidx_end, L, II):
    """compute results and write them to stream (rst_h = rst_ah + (rst_sum_sigma_h / (rst_sum_sigma + EPS)))

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_ah : simpy.store
        Input FIFO of result ah feature
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_in_rst_sum_sigma_h : simpy.store
        Input FIFO of result sum_sigma_h feature
    fifo_in_rst_sum_sigma : simpy.store
        Input FIFO of result sum_sigma feature
    fifo_out_rst_h : simpy.store
        Output FIFO of each 
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for i in range(2):
            if i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)

            if deg != 0:
                _ = yield fifo_in_rst_ah.get()
                _ = yield fifo_in_rst_sum_sigma_h.get()
                _ = yield fifo_in_rst_sum_sigma.get()
                yield fifo_out_rst_h.put(1)
            else:
                _ = yield fifo_in_rst_ah.get()
                yield fifo_out_rst_h.put(1)

    print("comp_rst ends at %d" % env.now)

def write_rst_h_mem(env, fifo_in_rst_h, nidx_begin, nidx_end, L, II, L_mem):
    """write node feature results to memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_h : simpy.store
        Input FIFO of result feature of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory write access
    """

    for n in range(nidx_begin, nidx_end):

        # yield env.timeout(L_mem) # memory latency

        for i in range(2):
            if i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)

            _ = yield fifo_in_rst_h.get()

    print("write_rst_h_mem ends at %d" % env.now)


if __name__ == "__main__":
    print("Starting GatedGCN HLPerf Model.")

    parser = argparse.ArgumentParser()
    add_option(parser)
    args = parser.parse_args()

    # read input graph
    format = "gnnhls"
    # input_fn = "./test/csr_indptr_trans.txt"
    # input_fn = "./test/csr_indptr_trans_v2.txt"

    # input_fn = "../GraphHLS/gcn_c/data_OGBG-MOLHIV/csr_indptr_trans.txt"
    # input_fn = "../GraphHLS/gcn_c/data_OGBN-PROTEINS/csr_indptr_trans.txt"

    # input_fn = "./data/OGBG-MOLHIV/csr_indptr_trans.txt"
    # input_fn = "./data/OGBN-PROTEINS/csr_indptr_trans.txt"
    # input_fn = "./data/OGBN-ARXIV/csr_indptr_trans.txt"
    # input_fn = "./data/OGBG-MOLTOX21/csr_indptr_trans.txt"
    input_fn = gen_input_fn(args.ds_fn)
    print("Input dataset:", args.ds_fn, ":", input_fn)

    # hardware parameters (MHz)
    freq = 295.4
    print("Clock frequency (MHz): ", freq)

    # generate degree list of the input graph
    deg_list = load_graph_deg(input_fn, format)
    # deg_list = [2] * 5000 + [100] * 5000
    # deg_list = [2, 100] * 5000

    # calculate node number and edge number
    node_num = len(deg_list)
    edge_num = sum(deg_list)
    print("node #:", node_num, "edge #:", edge_num)

    nidx_begin = 0
    nidx_end = node_num

    ft_size = 32
    d = 8
    print("feature size:", ft_size, "d:", d)


    start = timeit.default_timer()

    env = simpy.Environment()

    # read nod_src from memory
    fifo_nod_src_list = [simpy.Store(env, capacity=10), simpy.Store(env, capacity=11), simpy.Store(env, capacity=14), simpy.Store(env, capacity=15), simpy.Store(env, capacity=16), simpy.Store(env, capacity=17), simpy.Store(env, capacity=18), simpy.Store(env, capacity=19), simpy.Store(env, capacity=20), simpy.Store(env, capacity=21), simpy.Store(env, capacity=22), simpy.Store(env, capacity=23), simpy.Store(env, capacity=24), simpy.Store(env, capacity=25), simpy.Store(env, capacity=26), simpy.Store(env, capacity=27)]
    # fifo_nod_src_list = [simpy.Store(env, capacity=x+10) for x in range(15)]
    # print("len(fifo_nod_src_list):", len(fifo_nod_src_list))
    L = 4
    II = 1
    L_mem = 64
    env.process(read_nod_src(env, fifo_nod_src_list, nidx_begin, nidx_end, L, II, L_mem, deg_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L_mem + L + II * (node_num - 1)
    print("Estimated exec time of read_nod_src alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read edge_src from memory
    fifo_tmp_src = simpy.Store(env, capacity=15)
    L = 75
    II = 1
    L_mem = 64
    env.process(read_edge_src(env, fifo_nod_src_list[0], fifo_tmp_src, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3 + L_mem) * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of read_edge_src alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # read ft_in for layer A and E from memory (target features)
    fifo_ft_in_ah = simpy.Store(env, capacity=20)
    fifo_ft_in_eh = simpy.Store(env, capacity=26)
    L = 74
    II = 1
    L_mem = 64
    env.process(read_ft_in_tar(env, fifo_nod_src_list[1], fifo_ft_in_ah, fifo_ft_in_eh, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3 + L_mem + L + II) * node_num
    print("Estimated exec time of read_ft_in_tar alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

   
    # execute linear A layer - phase 1 (linear projection: in_feats -> out_feats)
    fifo_rst_ah_p1 = simpy.Store(env, capacity=10)
    L_list = [2, 13, 2]
    II_list = [1, 1, 1]
    env.process(update_tar_ah(env, fifo_ft_in_ah, fifo_rst_ah_p1, nidx_begin, nidx_end, L_list, II_list, ft_size, d))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + II_list[0] + L_list[1] + (ft_size-1) * II_list[1] + L_list[2] + II_list[2]) * node_num
    print("Estimated exec time of update_tar_ah alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # execute linear A layer - phase 2 (linear projection: in_feats -> out_feats)
    fifo_rst_ah = simpy.Store(env, capacity=48)
    L_list = [2, 9, 2]
    II_list = [1, 2, 1]
    env.process(update_tar_ah_sum(env, fifo_rst_ah_p1, fifo_rst_ah, nidx_begin, nidx_end, L_list, II_list, ft_size, d))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + II_list[0] + L_list[1] + (ft_size * d / 8 - 1) * II_list[1] + L_list[2] + II_list[2]) * node_num
    print("Estimated exec time of update_tar_ah_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # execute linear E layer - phase 1 (linear projection: in_feats -> out_feats)
    fifo_rst_eh_p1 = simpy.Store(env, capacity=10)
    L_list = [2, 13, 2]
    II_list = [1, 1, 1]
    env.process(update_tar_eh(env, fifo_ft_in_eh, fifo_nod_src_list[2], fifo_rst_eh_p1, nidx_begin, nidx_end, L_list, II_list, ft_size, d))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3 + L_list[0] + II_list[0] + L_list[1] + (ft_size-1) * II_list[1] + L_list[2] + (d - 1) * II_list[2]) * node_num
    print("Estimated exec time of update_tar_eh alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # execute linear E layer - phase 2 (linear projection: in_feats -> out_feats)
    fifo_rst_eh = simpy.Store(env, capacity=36)
    L_list = [2, 9, 2]
    II_list = [1, 2, 1]
    env.process(update_tar_eh_sum(env, fifo_rst_eh_p1, fifo_nod_src_list[3], fifo_rst_eh, nidx_begin, nidx_end, L_list, II_list, ft_size, d))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3 + L_list[0] * (d - 1) * II_list[0] + L_list[1] + (int(ft_size * d / 8) - 1) * II_list[1] + L_list[2] + II_list[2]) * node_num
    print("Estimated exec time of update_tar_eh_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # read aggregated features from memory
    fifo_ft_in_dh = simpy.Store(env, capacity=20)
    fifo_ft_in_bh = simpy.Store(env, capacity=24)
    L = 74
    II = 1
    L_mem = 64
    env.process(read_feat_in_agg(env, fifo_nod_src_list[4], fifo_tmp_src, fifo_ft_in_dh, fifo_ft_in_bh, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + (3) * edge_num + (L + II) * edge_num
    print("Estimated exec time of read_feat_in_agg alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # execute linear D layer - phase 1 (linear projection: in_feats -> out_feats)
    fifo_rst_dh_p1 = simpy.Store(env, capacity=10)
    L_list = [2, 13, 2]
    II_list = [1, 1, 1]
    env.process(update_dh(env, fifo_ft_in_dh, fifo_nod_src_list[5], fifo_rst_dh_p1, nidx_begin, nidx_end, L_list, II_list, ft_size, d))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + (L_list[0] + II_list[0] + L_list[1] + (ft_size-1) * II_list[1] + L_list[2] + (d - 1) * II_list[2]) * edge_num
    print("Estimated exec time of update_dh alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # execute linear D layer - phase 2 (linear projection: in_feats -> out_feats)
    fifo_rst_dh = simpy.Store(env, capacity=30)
    L_list = [2, 9, 2]
    II_list = [1, 2, 1]
    env.process(update_dh_sum(env, fifo_rst_dh_p1, fifo_nod_src_list[6], fifo_rst_dh, nidx_begin, nidx_end, L_list, II_list, ft_size, d))
    exec_time_alone = 3 * node_num + (L_list[0] + (d - 1) * II_list[0] + L_list[1] + (int(ft_size * d / 8) - 1) * II_list[1] + L_list[2] + II_list[2]) * edge_num
    print("Estimated exec time of update_dh_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # execute linear B layer - phase 1 (linear projection: in_feats -> out_feats)
    fifo_rst_bh_p1 = simpy.Store(env, capacity=5)
    L_list = [2, 13, 2]
    II_list = [1, 1, 1]
    env.process(update_bh(env, fifo_ft_in_bh, fifo_nod_src_list[7], fifo_rst_bh_p1, nidx_begin, nidx_end, L_list, II_list, ft_size, d))
    exec_time_alone = 3 * node_num + (L_list[0] + II_list[0] + L_list[1] + (ft_size-1) * II_list[1] + L_list[2] + (d - 1) * II_list[2]) * edge_num
    print("Estimated exec time of update_bh alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # execute linear B layer - phase 2 (linear projection: in_feats -> out_feats)
    fifo_rst_bh = simpy.Store(env, capacity=32)
    L_list = [2, 9, 2]
    II_list = [1, 2, 1]
    env.process(update_bh_sum(env, fifo_rst_bh_p1, fifo_nod_src_list[8], fifo_rst_bh, nidx_begin, nidx_end, L_list, II_list, ft_size, d))
    exec_time_alone = 3 * node_num + (L_list[0] + (d - 1) * II_list[0] + L_list[1] + (int(ft_size * d / 8) - 1) * II_list[1] + L_list[2] + II_list[2]) * edge_num
    print("Estimated exec time of update_bh_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read e_idx form memeory
    fifo_e_idx_list = [simpy.Store(env, capacity=2), simpy.Store(env, capacity=14)]
    L = 75
    II = 1
    L_mem = 64
    env.process(read_e_idx(env, fifo_nod_src_list[9], fifo_e_idx_list, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of read_e_idx alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read aggregated edge features from memory
    fifo_ft_in_ce = simpy.Store(env, capacity=20)
    L = 74
    II = 1
    L_mem = 64
    env.process(read_ft_in_edge_agg(env, fifo_nod_src_list[10], fifo_e_idx_list[0], fifo_ft_in_ce, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + (L + II) * edge_num
    print("Estimated exec time of read_ft_in_edge_agg alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # execute linear C layer: edge - phase 1 (linear projection: in_feats -> out_feats)
    fifo_rst_ce_p1 = simpy.Store(env, capacity=10)
    L_list = [2, 13, 2]
    II_list = [1, 1, 1]
    env.process(update_ce(env, fifo_ft_in_ce, fifo_nod_src_list[11], fifo_rst_ce_p1, nidx_begin, nidx_end, L_list, II_list, ft_size, d))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + (L_list[0] + II_list[0] + L_list[1] + (ft_size-1) * II_list[1] + L_list[2] + (d - 1) * II_list[2]) * edge_num
    print("Estimated exec time of update_ce alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # execute linear C layer: edge - phase 2 (linear projection: in_feats -> out_feats)
    fifo_rst_ce = simpy.Store(env, capacity=20)
    L_list = [2, 9, 2]
    II_list = [1, 2, 1]
    env.process(update_ce_sum(env, fifo_rst_ce_p1, fifo_nod_src_list[12], fifo_rst_ce, nidx_begin, nidx_end, L_list, II_list, ft_size, d))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + (L_list[0] + (d - 1) * II_list[0] + L_list[1] + (int(ft_size * d / 8) - 1) * II_list[1] + L_list[2] + II_list[2]) * edge_num
    print("Estimated exec time of update_ce_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # calculate rst_e, sum_sigma and sum_sigma_h
    fifo_rst_e = simpy.Store(env, capacity=20)
    fifo_rst_sum_sigma = simpy.Store(env, capacity=21)
    fifo_rst_sum_sigma_h = simpy.Store(env, capacity=21)
    L_list = [2, 65, 2]
    II_list = [1, 3, 1]
    env.process(comp_rst_e_sum_sigma(env, fifo_nod_src_list[13], fifo_rst_bh, fifo_rst_dh, fifo_rst_eh, fifo_rst_ce, fifo_rst_e, fifo_rst_sum_sigma, fifo_rst_sum_sigma_h, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3 + L_list[0] + II_list[0]) * node_num + L_list[1] * node_num + II_list[1] * (2 * edge_num - node_num) + (L_list[2] + II_list[2]) * node_num
    print("Estimated exec time of comp_rst_e_sum_sigma alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # write rst_e to memory
    L = 73
    II = 1
    L_mem = 64
    env.process(write_rst_e_mem(env, fifo_e_idx_list[1], fifo_nod_src_list[14], fifo_rst_e, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + (1 + L + II) * edge_num
    print("Estimated exec time of write_rst_e_mem alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # compute results and write them to stream (rst_h = rst_ah + (rst_sum_sigma_h / (rst_sum_sigma + EPS)))
    fifo_rst_h = simpy.Store(env, capacity=20)
    L = 29
    II = 1
    env.process(comp_rst(env, fifo_rst_ah, fifo_nod_src_list[15], fifo_rst_sum_sigma_h, fifo_rst_sum_sigma, fifo_rst_h, nidx_begin, nidx_end, L, II))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + (L + II) * node_num
    print("Estimated exec time of comp_rst alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # write node feature results to memory
    L = 74
    II = 1
    L_mem = 64
    env.process(write_rst_h_mem(env, fifo_rst_h, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L + II) * node_num
    print("Estimated exec time of write_rst_h_mem alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # run the simulation
    env.run()



    # measure the actural running time (seconds)
    end = timeit.default_timer()
    print("Actual running time: %s Seconds"%(end-start))

    # measure the simulation resultes (cycles)
    print("Simulation ends at %d" % env.now)
    print("Calculated simulation time (s): ", (env.now * 1000 / freq * 1e-9))