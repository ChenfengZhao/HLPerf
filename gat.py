import simpy
from utils import *
import timeit
import argparse

# ***** Kernel 1 *****
def read_feat_in(env, fifo_out_ft_in_tar, nidx_begin, nidx_end, L, II, L_mem):
    """read input feature from memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_out_ft_in_tar : simpy.store
        Output FIFO of input feature of each node
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

    print("Simulation starts at %d" % env.now)

    for n in range(nidx_begin, nidx_end):

        # yield env.timeout(L_mem)

        for i in range(2):
            if n == nidx_begin and i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            yield fifo_out_ft_in_tar.put(1)

    print("read_feat_in ends at %d" % env.now)

def update_tar_fc(env, fifo_in_ft_in_tar, fifo_out_rst_tar_p1, nidx_begin, nidx_end, L_list, II_list):
    """compute fc layer (multiplication) and write the calculated fc feature (ft_fc_array) to a stream channel (ft_fc_chan) (phase 1)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_in_tar : simpy.store
        Input FIFO of input feature of each node
    fifo_out_rst_tar_p1 : simpy.store
        Output FIFO of result feature (phase 1) of each node
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

        for i in range(2):
            _ = yield fifo_in_ft_in_tar.get()
            if i == 0:
                yield env.timeout(L_list[0])
            else:
                yield env.timeout(II_list[0])
        
        yield env.timeout(L_list[1] + (64-1) * II_list[1])

        for kd in range(16):
            if kd == 0:
                yield env.timeout(L_list[2])
            else:
                yield env.timeout(II_list[2])
            
            yield fifo_out_rst_tar_p1.put(1)

    print("update_tar_fc ends at %d" % env.now)

def update_tar_sum(env, fifo_in_rst_tar_p1, fifo_out_rst_fc, nidx_begin, nidx_end, L_list, II_list):
    """compute fc layer (multiplication) and write the calculated fc feature (ft_fc_array) to a stream channel (ft_fc_chan) (phase 2)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_tar_p1 : simpy.store
        Input FIFO of result feature (phase 1) of each node
    fifo_out_rst_fc : simpy.store
        Output FIFO of FC result feature of each node
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
        for kd in range(16):
            _ = yield fifo_in_rst_tar_p1.get()
            if kd == 0:
                yield env.timeout(L_list[0])
            else:
                yield env.timeout(II_list[0])
        
        yield env.timeout(L_list[1] + (128-1)*II_list[1])

        for i in range(8):
            if i == 0:
                yield env.timeout(L_list[2])
            else:
                yield env.timeout(II_list[2])
            
            yield fifo_out_rst_fc.put(1)

    print("update_tar_sum ends at %d" % env.now)

def split1(env, fifo_in_rst_fc, fifo_out_rst_fc_branch_list, nidx_begin, nidx_end, L, II):
    """read stream channel ft_fc_chan to 3 other stream channels

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_fc : simpy.store
        Input FIFO of FC result feature of each node
    fifo_out_rst_fc_branch_list : list of simpy.store
        Output FIFO list of FC result feature copys of each node
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
        for i in range(8):
            _ = yield fifo_in_rst_fc.get()

            if n == nidx_begin and i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)

            for j in range(len(fifo_out_rst_fc_branch_list)):

                yield fifo_out_rst_fc_branch_list[j].put(1)
    
    print("split1 ends at %d" % env.now)

def write_feat_fc(env, fifo_in_rst_fc_branch, nidx_begin, nidx_end, L, II, L_mem):
    """write fc results to memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_fc_branch : simpy.store
        Input FIFO of FC result feature of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    L_mem : int
        The latency of memory write access
    """

    for n in range(nidx_begin, nidx_end):
        for i in range(8):
            _ = yield fifo_in_rst_fc_branch.get()
            
            if n == nidx_begin and i == 0:
                yield env.timeout(L+L_mem)
            else:
                yield env.timeout(II)
    
    print("write_feat_fc ends at %d" % env.now)

def comp_el(env, fifo_in_rst_fc_branch, fifo_out_rst_el, nidx_begin, nidx_end, L, II):
    """compute el

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_fc_branch : simpy.store
        Input FIFO of FC result feature of each node
    fifo_out_rst_el : simpy.store
        Output FIFO of result el feature of each node
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
        for i in range(8):
            _ = yield fifo_in_rst_fc_branch.get()

            if n == nidx_begin and i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            yield fifo_out_rst_el.put(1)

    print("comp_el ends at %d" % env.now)

def write_el_mem(env, fifo_in_rst_el, nidx_begin, nidx_end, L, II, L_mem):
    """write el from stream channel to memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_el : simpy.store
        Input FIFO of calcutaed el feature of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        The latency of memory write access
    """

    for n in range(nidx_begin, nidx_end):
        for i in range(8):
            _ = yield fifo_in_rst_el.get()

            if n == nidx_begin and i == 0:
                yield env.timeout(L+L_mem)
            else:
                yield env.timeout(II)
    
    print("write_el_mem ends at %d" % env.now)

def comp_er(env, fifo_in_rst_fc_branch, fifo_out_rst_er, nidx_begin, nidx_end, L, II):
    """compute er

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_fc_branch : simpy.store
        Input FIFO of FC result feature of each node
    fifo_out_rst_er : simpy.store
        Output FIFO of calculated er feature of each node
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
        for i in range(8):
            _ = yield fifo_in_rst_fc_branch.get()

            if n == nidx_begin and i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            yield fifo_out_rst_er.put(1)
    
    print("comp_er ends at %d" % env.now)

def write_er_mem(env, fifo_in_rst_er, nidx_begin, nidx_end, L, II, L_mem):
    """write er from stream channer to memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_er : simpy.store
        FIFO Input calculated er feature of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        The latency of memory write access
    """

    for n in range(nidx_begin, nidx_end):
        for i in range(8):
            _ = yield fifo_in_rst_er.get()

            if n == nidx_begin and i == 0:
                yield env.timeout(L+L_mem)
            else:
                yield env.timeout(II)
    
    print("write_er_mem ends at %d" % env.now)

# ***** Kernel 2 *****
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

def read_edge_src(env, fifo_in_nod_src, fifo_out_tmp_src_list, nidx_begin, nidx_end, L, II, L_mem):
    """read edge_src from memory
    
    Instead of putting actual edge index into the FIFO, we put fake values

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_tmp_src_list : List of simpy.store
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
            
            for i in range(len(fifo_out_tmp_src_list)):
                yield fifo_out_tmp_src_list[i].put(1)

    print("read_edge_src ends at %d" % env.now)

def read_mem_er(env, fifo_in_nod_src, fifo_out_er, nidx_begin, nidx_end, L, II, L_mem):
    """read er_mat from memory to local buffer (er_array)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_er : simpy.store
        Output FIFO of er feature of each node
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
        yield env.timeout(4)

        if deg != 0:

            yield env.timeout(L_mem)

            for i in range(8):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)
                
                yield fifo_out_er.put(1)

    print("read_mem_er ends at %d" % env.now)  

def read_mem_el(env, fifo_in_nod_src, fifo_in_tmp_src, fifo_out_el, nidx_begin, nidx_end, L, II, L_mem):
    """read el from memory to stream channel

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_in_tmp_src : simpy.store
        Input FIFO of neighbor indices of each node
    fifo_out_el : simpy.store
        Output FIFO of el feature of each neighbor
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
            yield env.timeout(1)
            _ = yield fifo_in_tmp_src.get()

            yield env.timeout(L_mem) # memory latency

            for i in range(8):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)
                
                yield fifo_out_el.put(1)

    print("read_mem_el ends at %d" % env.now)

def read_edge_src2(env, fifo_in_nod_src, fifo_out_tmp_src, nidx_begin, nidx_end, L, II, L_mem):
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

    print("read_edge_src2 ends at %d" % env.now)

# read_mem_er_2 is the same as read_mem_er
def read_mem_er2(env, fifo_in_nod_src, fifo_out_er, nidx_begin, nidx_end, L, II, L_mem):
    """read er_mat from memory to local buffer (er_array)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_er : simpy.store
        Output FIFO of er feature of each node
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
        yield env.timeout(4)

        if deg != 0:

            yield env.timeout(L_mem)

            for i in range(8):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)
                
                yield fifo_out_er.put(1)

    print("read_mem_er2 ends at %d" % env.now)  

# read_mem_el_2 is the same as read_mem_el
def read_mem_el2(env, fifo_in_nod_src, fifo_in_tmp_src, fifo_out_el, nidx_begin, nidx_end, L, II, L_mem):
    """read el from memory to stream channel

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_in_tmp_src : simpy.store
        Input FIFO of neighbor indices of each node
    fifo_out_el : simpy.store
        Output FIFO of el feature of each neighbor
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
            yield env.timeout(1)
            _ = yield fifo_in_tmp_src.get()

            yield env.timeout(L_mem) # memory latency

            for i in range(8):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)
                
                yield fifo_out_el.put(1)

    print("read_mem_el2 ends at %d" % env.now)

def comp_e(env, fifo_in_el, fifo_in_er, fifo_in_nod_src, fifo_out_e, nidx_begin, nidx_end, L_list, II_list):
    """compute e = leaky_relu(el + er) for all heads and write computed e to stream channel (e_array)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_el : simpy.store
        Input FIFO of el feature of each neighbor
    fifo_in_er : simpy.store
        Input FIFO of er feature of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_e : simpy.store
        Output feture of calculated e feature of each edge/neighbor
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
        yield env.timeout(4)

        if deg != 0:

            for i in range(8):
                _ = yield fifo_in_er.get()

                if i == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            for e in range(deg):
                for i in range(8):
                    _ = yield fifo_in_el.get()

                    if e == 0 and i == 0:
                        yield env.timeout(L_list[1])
                    else:
                        yield env.timeout(II_list[1])
                    
                    yield fifo_out_e.put(1)

    print("comp_e ends at %d" % env.now)

def comp_e_sum(env, fifo_in_el, fifo_in_er, fifo_in_nod_src, fifo_out_ek_sum, nidx_begin, nidx_end, L_list, II_list):
    """compute e = leaky_relu(el + er) for all heads, aggregate computed e for each node, and write computed e to stream channel (e_array)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_el : simpy.store
        Input FIFO of el feature of each neighbor
    fifo_in_er : simpy.store
        Input FIFO of er feature of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_ek_sum : simpy.store
        Output FIFO of aggregated ek_sum feature of each node
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
        yield env.timeout(4)

        if deg != 0:

            for i in range(8):
                _ = yield fifo_in_er.get()

                if i == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])

            for e in range(deg):
                for i in range(8):
                    _ = yield fifo_in_el.get()

                    if e == 0 and i == 0:
                        yield env.timeout(L_list[1])
                    else:
                        yield env.timeout(II_list[1])
            
            for i in range(8):
                if i == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
            
                yield fifo_out_ek_sum.put(1)

    print("comp_e_sum ends at %d" % env.now)

def comp_softmax(env, fifo_in_e, fifo_in_ek_sum, fifo_in_nod_src, fifo_out_a, nidx_begin, nidx_end, L_list, II_list):
    """compute softmax

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_e : simpy.store
        Input FIFO of calculated e feature of each edge/neighbor
    fifo_in_ek_sum : simpy.store
        Input FIFO of aggregated ek_sum feature of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_a : simpy.store
        Output FIFO of softmax results of each edge
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
        yield env.timeout(4)

        if deg != 0:

            for i in range(8):
                _ = yield fifo_in_ek_sum.get()
                if i == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            for e in range(deg):
                for i in range(8):
                    _ = yield fifo_in_e.get()

                    if e == 0 and i == 0:
                        yield env.timeout(L_list[1])
                    else:
                        yield env.timeout(II_list[1])
                    
                    yield fifo_out_a.put(1)

    print("comp_softmax ends at %d" % env.now)


def read_mem_fc_nbr(env, fifo_in_nod_src, fifo_in_tmp_src, fifo_out_ft_fc_nbr, nidx_begin, nidx_end, L, II, L_mem):
    """read neighbor input feature from memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_in_tmp_src : simpy.store
        Input FIFO of neighbor indices of each node
    fifo_out_ft_fc_nbr : simpy.store
        Output FIFO of features of each neighbor/edge
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
            yield env.timeout(1)
            _ = yield fifo_in_tmp_src.get()

            # yield env.timeout(L_mem) # memory latency

            for i in range(8):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)
                
                yield fifo_out_ft_fc_nbr.put(1)

    print("read_mem_fc_nbr ends at %d" % env.now)

def comp_rst(env, fifo_in_ft_fc_nbr, fifo_in_a, fifo_in_nod_src, fifo_out_rst, nidx_begin, nidx_end, L_list, II_list):
    """compute final results: message passing (aggregate (ft_fc * a) and sum)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_fc_nbr : simpy.store
        Input FIFO of features of each neighbor/edge
    fifo_in_a : simpy.store
        Input FIFO of softmax results of each edge
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst : simpy.store
        Output FIFO of final results feature of each node
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
        yield env.timeout(4)

        if deg != 0:

            for e in range(deg):
                for i in range(8):
                    _ = yield fifo_in_a.get()
                    _ = yield fifo_in_ft_fc_nbr.get()

                    if e == 0 and i == 0:
                        yield env.timeout(L_list[0])
                    else:
                        yield env.timeout(II_list[0])
            
            for i in range(8):
                if i == 0:
                    yield env.timeout(L_list[1])
                else:
                    yield env.timeout(II_list[1])
                
                yield fifo_out_rst.put(1)

    print("comp_rst ends at %d" % env.now)

def wirte_mem_rst(env, fifo_in_rst, fifo_in_nod_src, nidx_begin, nidx_end, L, II, L_mem):
    """write final results to memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst : simpy.store
        Input FIFO of final results feature of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node
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
        yield env.timeout(4)

        if deg != 0:

            yield env.timeout(L_mem) # memory latency

            for i in range(8):
                _ = yield fifo_in_rst.get()

                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)

    print("wirte_mem_rst ends at %d" % env.now)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    add_option(parser)
    args = parser.parse_args()

    # read input graph
    format = "gnnhls"
    # input_fn = "./test/csr_indptr_trans.txt"
    # input_fn = "./test/csr_indptr_trans_v2.txt"

    # input_fn = "../GraphHLS/gcn_c/data_OGBG-MOLHIV/csr_indptr_trans.txt"
    # input_fn = "../GraphHLS/gcn_c/data_OGBN-PROTEINS/csr_indptr_trans.txt"

    # input_fn = "./data/OGBG-MOLTOX21/csr_indptr_trans.txt"
    # input_fn = "./data/OGBG-MOLHIV/csr_indptr_trans.txt"
    # input_fn = "./data/OGBN-ARXIV/csr_indptr_trans.txt"
    # input_fn = "./data/OGBN-PROTEINS/csr_indptr_trans.txt"

    # input_fn = input_fn_dict[args.ds_fn]
    input_fn = gen_input_fn(args.ds_fn)
    print("Input dataset:", args.ds_fn, ":", input_fn)

    # hardware parameters (MHz)
    freq = 225.9
    print("CLK Frequency:", freq, "MHz")

    # generate degree list of the input graph
    deg_list = load_graph_deg(input_fn, format)
    # deg_list = [2] * 5000 + [50] * 5000
    # deg_list = [2, 50] * 5000

    # calculate node number and edge number
    node_num = len(deg_list)
    edge_num = sum(deg_list)
    print("node #:", node_num, "edge #:", edge_num)

    nidx_begin = 0
    nidx_end = node_num

    print("***** Kernel 1 *****")
    start = timeit.default_timer()

    env1 = simpy.Environment()

    # read input feature from memory
    fifo_ft_in_tar = simpy.Store(env1, capacity=10)
    L = 77
    II = 2
    L_mem = 64
    env1.process(read_feat_in(env1, fifo_ft_in_tar, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + (II) * node_num
    print("Estimated exec time of read_feat_in alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # compute fc layer (multiplication) and write the calculated fc feature (ft_fc_array) to a stream channel (ft_fc_chan) (phase 1)
    fifo_rst_tar_p1 = simpy.Store(env1, capacity=16)
    L_list = [2, 13, 2]
    II_list = [1, 1, 1]
    env1.process(update_tar_fc(env1, fifo_ft_in_tar, fifo_rst_tar_p1, nidx_begin, nidx_end, L_list, II_list))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + II_list[0]) * node_num + (L_list[1] + (64-1) * II_list[1]) * node_num + (L_list[2] + (16-1) * II_list[2]) * node_num
    print("Estimated exec time of update_tar_fc alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # compute fc layer (multiplication) and write the calculated fc feature (ft_fc_array) to a stream channel (ft_fc_chan) (phase 2)
    fifo_rst_fc = simpy.Store(env1, capacity=10)
    L_list = [2, 9, 2]
    II_list = [1, 1, 1]
    env1.process(update_tar_sum(env1, fifo_rst_tar_p1, fifo_rst_fc, nidx_begin, nidx_end, L_list, II_list))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + (16-1) * II_list[0]) * node_num + (L_list[1] + (128-1)*II_list[1]) * node_num + (L_list[2] + (8-1) * II_list[2]) * node_num
    print("Estimated exec time of update_tar_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read stream channel ft_fc_chan to 3 other stream channels
    fifo_rst_fc_branch_list = [simpy.Store(env1, capacity=10), simpy.Store(env1, capacity=11), simpy.Store(env1, capacity=12)]
    L = 3
    II = 1
    env1.process(split1(env1, fifo_rst_fc, fifo_rst_fc_branch_list, nidx_begin, nidx_end, L, II))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + II * (node_num * 8 -1)
    print("Estimated exec time of split1 alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # write fc results to memory
    L = 4
    II = 1
    L_mem = 64
    env1.process(write_feat_fc(env1, fifo_rst_fc_branch_list[0], nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + L_mem + II * (node_num * 8 - 1)
    print("Estimated exec time of write_feat_fc alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # compute el
    fifo_rst_el = simpy.Store(env1, capacity=16)
    L = 113
    II = 1
    env1.process(comp_el(env1, fifo_rst_fc_branch_list[1], fifo_rst_el, nidx_begin, nidx_end, L, II))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + II * (node_num * 8 - 1)
    print("Estimated exec time of comp_el alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # write el from stream channel to memory
    L = 4
    II = 1
    L_mem = 64
    env1.process(write_el_mem(env1, fifo_rst_el, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + L_mem + II * (node_num * 8 - 1)
    print("Estimated exec time of write_el_mem alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # compute er
    fifo_rst_er = simpy.Store(env1, capacity=16)
    L = 113
    II = 1
    env1.process(comp_er(env1, fifo_rst_fc_branch_list[2], fifo_rst_er, nidx_begin, nidx_end, L, II))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + II * (node_num * 8 - 1)
    print("Estimated exec time of comp_er alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # write er from stream channer to memory
    L = 4
    II = 1
    L_mem = 64
    env1.process(write_er_mem(env1, fifo_rst_er, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + L_mem + II * (node_num * 8 - 1)
    print("Estimated exec time of write_er_mem alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # run the simulation
    env1.run()

    # measure the actural running time (seconds)
    end = timeit.default_timer()
    print("Actual running time: %s Seconds"%(end-start))

    # measure the simulation resultes (cycles)
    print("Simulation ends at %d" % env1.now)
    print("Calculated simulation time (s): ", (env1.now * 1000 / freq * 1e-9))



    print("***** Kernel 2 *****")

    nidx_begin = 0
    nidx_end = node_num

    start = timeit.default_timer()

    env2 = simpy.Environment()

    # read nod_src from memory
    # fifo_nod_src_list = [simpy.Store(env2, capacity=10), simpy.Store(env2, capacity=11), simpy.Store(env2, capacity=12), simpy.Store(env2, capacity=13), simpy.Store(env2, capacity=14), simpy.Store(env2, capacity=15), simpy.Store(env2, capacity=16), simpy.Store(env2, capacity=17), simpy.Store(env2, capacity=18), simpy.Store(env2, capacity=19), simpy.Store(env2, capacity=20), simpy.Store(env2, capacity=21)]
    fifo_nod_src_list = [simpy.Store(env2, capacity=(10+i)) for i in range(12)]
    # fifo_nod_src_list = [simpy.Store(env2, capacity=10), simpy.Store(env2, capacity=11)]
    L = 4
    II = 1
    L_mem = 64
    env2.process(read_nod_src(env2, fifo_nod_src_list, nidx_begin, nidx_end, L, II, L_mem, deg_list))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L_mem + L + II * (node_num - 1)
    print("Estimated exec time of read_nod_src alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # read edge_src from memory
    fifo_tmp_src_list = [simpy.Store(env2, capacity=10), simpy.Store(env2, capacity=17)]
    # fifo_tmp_src_list = [simpy.Store(env2, capacity=10)]
    L = 75
    II = 1
    L_mem = 64
    env2.process(read_edge_src(env2, fifo_nod_src_list[0], fifo_tmp_src_list, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3 + L_mem) * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of read_edge_src alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # read er_mat from memory to local buffer (er_array)
    fifo_er = simpy.Store(env2, capacity=80)
    L = 2
    II = 1
    L_mem = 64
    env2.process(read_mem_er(env2, fifo_nod_src_list[1], fifo_er, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + L_mem * node_num + L * node_num + II * ((8-1) * node_num)
    print("Estimated exec time of read_mem_er alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # read el from memory to stream channel
    fifo_el = simpy.Store(env2, capacity=80)
    L = 2
    II = 1
    L_mem = 64
    env2.process(read_mem_el(env2, fifo_nod_src_list[2], fifo_tmp_src_list[0], fifo_el, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + 1 * edge_num + L_mem * edge_num + L * edge_num + II * ((8-1) * edge_num)
    print("Estimated exec time of read_mem_el alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    
    # compute e = leaky_relu(el + er) for all heads and write computed e to stream channel (e_array)
    fifo_e = simpy.Store(env2, capacity=80)
    L_list = [2, 30]
    II_list = [1, 1]
    env2.process(comp_e(env2, fifo_el, fifo_er, fifo_nod_src_list[3], fifo_e, nidx_begin, nidx_end, L_list, II_list))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + L_list[0] * node_num + II_list[0] * (8 - 1) * node_num + L_list[1] * node_num + II_list[1] * (8 * edge_num - node_num)
    print("Estimated exec time of comp_e alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))


    # read edge_src2 from memory for sum
    fifo_tmp_src2 = simpy.Store(env2, capacity=11)
    L = 75
    II = 1
    L_mem = 64
    env2.process(read_edge_src2(env2, fifo_nod_src_list[4], fifo_tmp_src2, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of read_edge_src2 alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read er_mat2 from memory to local buffer for sum (read_mem_er_2)
    fifo_er2 = simpy.Store(env2, capacity=80)
    L = 2
    II = 1
    L_mem = 64
    env2.process(read_mem_er2(env2, fifo_nod_src_list[5], fifo_er2, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + L_mem * node_num + L * node_num + II * ((8-1) * node_num)
    print("Estimated exec time of read_mem_er2 alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read el_mat2 from memory to local buffer for sum
    fifo_el2 = simpy.Store(env2, capacity=80)
    L = 2
    II = 1
    L_mem = 64
    env2.process(read_mem_el2(env2, fifo_nod_src_list[6], fifo_tmp_src2, fifo_el2, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + 1 * edge_num + L_mem * edge_num + L * edge_num + II * ((8-1) * edge_num)
    print("Estimated exec time of read_mem_el2 alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # compute e_sum
    fifo_ek_sum = simpy.Store(env2, capacity=80)
    L_list = [2, 35, 2]
    II_list = [1, 1, 1]
    env2.process(comp_e_sum(env2, fifo_el2, fifo_er2, fifo_nod_src_list[7], fifo_ek_sum, nidx_begin, nidx_end, L_list, II_list))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + L_list[0] * node_num + II_list[0] * ((8 - 1) * node_num) + L_list[1] * node_num + II_list[1] * (edge_num * 8 - node_num) + L_list[2] * node_num + II_list[2] * ((8 - 1) * node_num)
    print("Estimated exec time of comp_e_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # compute softmax
    fifo_a = simpy.Store(env2, capacity=80)
    L_list = [2, 16]
    II = [1, 1]
    env2.process(comp_softmax(env2, fifo_e, fifo_ek_sum, fifo_nod_src_list[8], fifo_a, nidx_begin, nidx_end, L_list, II_list))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + L_list[0] * node_num + II_list[0] * ((8-1) * node_num) + L_list[1] * node_num + II_list[1] * (edge_num * 8 - node_num)
    print("Estimated exec time of comp_softmax alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read neighbor input feature from memory
    fifo_ft_fc_nbr = simpy.Store(env2, capacity=80)
    L = 74
    II = 1
    L_mem = 64
    env2.process(read_mem_fc_nbr(env2, fifo_nod_src_list[9], fifo_tmp_src_list[1], fifo_ft_fc_nbr, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + (1 + L_mem) * edge_num + L * edge_num + II * ((8-1) * edge_num)
    print("Estimated exec time of read_mem_fc_nbr alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # compute final results: message passing (aggregate (ft_fc * a) and sum)
    fifo_rst = simpy.Store(env2, capacity=80)
    L_list = [14, 2]
    II_list = [1, 1]
    env2.process(comp_rst(env2, fifo_ft_fc_nbr, fifo_a, fifo_nod_src_list[10], fifo_rst, nidx_begin, nidx_end, L_list, II_list))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + L_list[0] * node_num + II_list[0] * (edge_num * 8 - node_num) + L_list[1] * node_num + II_list[1] * ((8-1) * node_num)
    print("Estimated exec time of comp_rst alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # write final results to memory
    L = 25
    II = 1
    L_mem = 64
    env2.process(wirte_mem_rst(env2, fifo_rst, fifo_nod_src_list[11], nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + L_mem * node_num + L * node_num + II * (8-1) * node_num
    print("Estimated exec time of wirte_mem_rst alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # run the simulation
    env2.run()

    # measure the actural running time (seconds)
    end = timeit.default_timer()
    print("Actual running time: %s Seconds"%(end-start))

    # measure the simulation resultes (cycles)
    print("Simulation ends at %d" % env2.now)
    print("Calculated simulation time (s): ", (env2.now * 1000 / freq * 1e-9))






