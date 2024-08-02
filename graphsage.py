import simpy
from utils import *
import timeit
import argparse

def read_nod_src(env, fifo_out_nod_src, nidx_begin, nidx_end, L, II, L_mem, deg_list):
    """read nod_src from memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_out_nod_src : simpy.store
        Output FIFO of neighbor index range of each node
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
        The list of node (in-comming) degree.
    """
    print("Simulation starts at %d" % env.now)

    for n in range(nidx_begin, nidx_end):

        if n == nidx_begin:
            yield env.timeout(L+L_mem)
        else:
            yield env.timeout(II)
        
        yield fifo_out_nod_src.put(deg_list[n])
    
    print("read_nod_src ends at %d" % env.now)

def split_nod_stream(env, fifo_in_nod_src, fifo_out_nod_src_list, nidx_begin, nidx_end, L, II):
    """split the nod_src_stream

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
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
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()

        if n == nidx_begin:
            yield env.timeout(L)
        else:
            yield env.timeout(II)

        for fifo_idx in range(len(fifo_out_nod_src_list)):
            yield fifo_out_nod_src_list[fifo_idx].put(deg)

    print("split_nod_stream ends at %d" % env.now)

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

def read_feat_in_agg(env, fifo_in_tmp_src, fifo_in_nod_src, fifo_out_ft_in_agg, nidx_begin, nidx_end, L, II, L_mem):
    """read input features to be aggregated from memory

    Instead of putting acutal input feature of nbrs, we put feke values

    Parameters  
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_tmp_src :simpy.store
        Input FIFO of neighbor indices of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_ft_in_agg : simpy.store
        Output FIFO of input feature to be aggretated
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

        for e in range(deg):
            yield env.timeout(1)
            _ = yield fifo_in_tmp_src.get()

            # yield env.timeout(L_mem) # memory latency
            # yield env.timeout(15) # memory latency

            for i in range(2):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)

                yield fifo_out_ft_in_agg.put(1)

    print("read_feat_in_agg ends at %d" % env.now)

def agg_feat_in(env, fifo_in_ft_in_agg, fifo_in_nod_src, fifo_out_ft_h_agg, nidx_begin, nidx_end, L_list, II_list):
    """aggregate features

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_in_agg : simpy.store
        Input FIFO of input feature to be aggretated
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_ft_h_agg : simpy.store
        Output FIFO of aggretaed fature
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency of the loop
    II_list : list of int
        The II of the loop
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        if deg != 0:
            for e in range(deg):

                for i in range(2):
                    if e == 0 and i == 0:
                        yield env.timeout(L_list[0])
                    else:
                        yield env.timeout(II_list[0])
                    _ = yield fifo_in_ft_in_agg.get()
            
            yield env.timeout(1)

            for i in range(2):
                if i == 0:
                    yield env.timeout(L_list[1])
                else:
                    yield env.timeout(II_list[1])
            
                yield fifo_out_ft_h_agg.put(1)
    
    print("agg_feat_in ends at %d" % env.now)



def update_agg(env, fifo_in_ft_h_agg, fifo_in_nod_src, fifo_out_rst_agg_p1, nidx_begin, nidx_end, L_list, II_list):
    """update/apply phase for the aggrated results (phase 1)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_h_agg : simpy.store
        Input FIFO of aggregated feature of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_agg_p1 : simpy.store
        Out FIFO of result feature (phase 1) of each node
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
        # print("n in update_agg:", n)
        deg = yield fifo_in_nod_src.get()
        # print("deg in update_agg:", deg)
        yield env.timeout(4)

        if deg != 0:
            for i in range(2):
                _ = yield fifo_in_ft_h_agg.get()
                if i == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
                
            # for k in range(64):
            #     if k == 0:
            #         yield env.timeout(L_list[1])
            #     else:
            #         yield env.timeout(II_list[1])
            
            yield env.timeout(L_list[1] + (64-1) * II_list[1])
            
            for kd in range(16):
                if kd == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
                
                yield fifo_out_rst_agg_p1.put(1)
        # else:
        #     yield env.timeout(160)

    print("update_agg ends at %d" % env.now)

def update_agg_sum(env, fifo_in_rst_agg_p1, fifo_in_nod_src, fifo_out_rst_agg, nidx_begin, nidx_end, L_list, II_list):
    """update/apply phase for the aggrated results (phase 2)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_agg_p1 : simpy.store
        Input FIFO of result feature (phase 1) of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_agg : simpy.store
        Output FIFO of result updated feature (phase 2) of each node
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
            for kd in range(16):
                _ = yield fifo_in_rst_agg_p1.get()
                if kd == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            # for kd in range(128):
            #     if kd == 0:
            #         yield env.timeout(L_list[1])
            #     else:
            #         yield env.timeout(II_list[1])
            
            yield env.timeout(L_list[1] + (128-1)*II_list[1])
            
            for i in range(2):
                if i == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
                yield fifo_out_rst_agg.put(1)

    print("update_agg_sum ends at %d" % env.now)


def read_feat_in_tar(env, fifo_out_ft_in_tar, nidx_begin, nidx_end, L, II, L_mem):
    """read the target feature from memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_out_ft_in_tar : simpy.store
        Output fifo of input feature of the target node
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
        # yield env.timeout(L_mem)

        for i in range(2):
            if n == nidx_begin and i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            yield fifo_out_ft_in_tar.put(1)

    print("read_feat_in_tar ends at %d" % env.now)


def update_tar(env, fifo_in_ft_in_tar, fifo_out_rst_tar_p1, nidx_begin, nidx_end, L_list, II_list):
    """update/apply phase for the target feature (phase 1)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_in_tar : simpy.store
        Input FIFO of the feature of each target node
    fifo_out_rst_tar_p1 : simpy.store
        Output FIFO of result feature (phase 1) of each target node
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

    print("update_tar ends at %d" % env.now)

def update_tar_sum(env, fifo_in_rst_tar_p1, fifo_out_rst_tar, nidx_begin, nidx_end, L_list, II_list):
    """update/apply phase for the target feature (phase 2)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_tar_p1 : simpy.store
        Input FIFO of result feature (phase 1) of each target node
    fifo_out_rst_tar : simpy.store
        Output FIFO of result feature (phase 2) of each target node
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

        for i in range(2):
            if i == 0:
                yield env.timeout(L_list[2])
            else:
                yield env.timeout(II_list[2])
            yield fifo_out_rst_tar.put(1)
    
    print("update_tar_sum ends at %d" % env.now)

def concat_update_rst(env, fifo_in_rst_agg, fifo_in_rst_tar, fifo_in_nod_src, fifo_out_rst_cat, nidx_begin, nidx_end, L, II):
    """concat the update/apply results of aggregated features and the target feature

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_agg : simpy.store
        Input FIFO of result aggregated feature
    fifo_in_rst_tar : simpy.store
        Input FIFO of result target feature
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_cat : simpy.store
        Ouput FIFO of result cat feature
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
                _ = yield fifo_in_rst_agg.get()
                _ = yield fifo_in_rst_tar.get()
                yield fifo_out_rst_cat.put(1)
            else:
                _ = yield fifo_in_rst_tar.get()
                yield fifo_out_rst_cat.put(1)

    print("concat_update_rst ends at %d" % env.now)

def norm_reduce_sum(env, fifo_in_rst_cat, fifo_out_rst_cat, fifo_out_rst_norm_reduce_sum, nidx_begin, nidx_end, L, II):
    """normalization for each vertex (cal base + norm div)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_cat : simpy.store
        Input FIFO of result cat feature
    fifo_out_rst_cat : simpy.store
        Output FIFO of result cat feature
    fifo_out_rst_norm_reduce_sum : _type_
        Output FIFO of noralized result cat feature
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
        for i in range(2):
            _ = yield fifo_in_rst_cat.get()

            if n == nidx_begin and i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            yield fifo_out_rst_cat.put(1)
            yield fifo_out_rst_norm_reduce_sum.put(1)

    print("norm_reduce_sum ends at %d" % env.now)

def norm_acc_rst(env, fifo_in_rst_norm_reduce_sum, fifo_out_rst_norm_acc, nidx_begin, nidx_end, L, II):
    """accumulate norm reduced sum results of all FIFO segments

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_norm_reduce_sum : simpy.store
        Input FIFO of noralized result cat feature
    fifo_out_rst_norm_acc : simpy.store
        Output FIFO of accumulated norm reduced sum feature
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
        for i in range(2):
            _ = yield fifo_in_rst_norm_reduce_sum.get()
            if i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
        yield env.timeout(10)
        yield fifo_out_rst_norm_acc.put(1)

    print("norm_acc_rst ends at %d" % env.now)

def norm_div(env, fifo_in_rst_cat, fifo_in_rst_norm_acc, fifo_out_rst_norm, nidx_begin, nidx_end, L, II):
    """divide concat results by norm accumulate results

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_cat : simpy.store
        Input FIFO of result cat feature
    fifo_in_rst_norm_acc : simpy.store
        Input FIFO of accumulated norm reduced sum feature
    fifo_out_rst_norm : simpy.store
        Output FIFO of divid concat result features
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
        for i in range(2):

            if i == 0:
                _ = yield fifo_in_rst_norm_acc.get()
            
            _ = yield fifo_in_rst_cat.get()

            if n == nidx_begin and i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)

            yield fifo_out_rst_norm.put(1)

    print("norm_div ends at %d" % env.now)

def write_norm_rst_mem(env, fifo_in_rst_norm, nidx_begin, nidx_end, L_list, II_list, L_mem):
    """write norm results to memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_norm : simpy.store
        Input FIFO of norm result feature
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list
        The list of loop latency 
    II_list : list
        The list of loop II
    L_mem : int
        the latency of memory write access
    """
    for n in range(nidx_begin, nidx_end):
        for i in range(2):
            _ = yield fifo_in_rst_norm.get()
            
            if i == 0:
                yield env.timeout(L_list[0])
            else:
                yield env.timeout(II_list[0])
        
        yield env.timeout(L_mem)

        for i in range(2):
            if i == 0:
                yield env.timeout(L_list[1])
            else:
                yield env.timeout(II_list[1])
    
    print("write_norm_rst_mem ends at %d" % env.now)

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
    freq = 165.1
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

    start = timeit.default_timer()

    env = simpy.Environment()

    # read nod_src from memory
    fifo_nod_src = simpy.Store(env, capacity=10)
    L = 4
    II = 1
    L_mem = 64
    env.process(read_nod_src(env, fifo_nod_src, nidx_begin, nidx_end, L, II, L_mem, deg_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L_mem + L + II * (node_num - 1)
    print("Estimated exec time of read_nod_src alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # split the nod_src_stream
    fifo_nod_src_list = [simpy.Store(env, capacity=11), simpy.Store(env, capacity=12), simpy.Store(env, capacity=13), simpy.Store(env, capacity=14), simpy.Store(env, capacity=15), simpy.Store(env, capacity=17)]
    L = 3
    II = 1
    env.process(split_nod_stream(env, fifo_nod_src, fifo_nod_src_list, nidx_begin, nidx_end, L, II))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + II * (node_num - 1)
    print("Estimated exec time of split_nod_stream alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read edge_src from memory
    fifo_tmp_src = simpy.Store(env, capacity=10)
    L = 75
    II = 1
    L_mem = 64
    env.process(read_edge_src(env, fifo_nod_src_list[0], fifo_tmp_src, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3 + L_mem) * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of read_edge_src alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read input features to be aggregated from memory
    fifo_ft_in_agg = simpy.Store(env, capacity=10)
    L = 74
    II = 1
    L_mem = 64
    env.process(read_feat_in_agg(env, fifo_tmp_src, fifo_nod_src_list[1], fifo_ft_in_agg, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + (L_mem + 1) * edge_num + (L + II) * edge_num
    print("Estimated exec time of read_feat_in_agg alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # aggregate features
    fifo_ft_h_agg = simpy.Store(env, capacity=10)
    L_list = [10, 15]
    II_list = [1, 1]
    env.process(agg_feat_in(env, fifo_ft_in_agg, fifo_nod_src_list[2], fifo_ft_h_agg, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + L_list[0] * node_num + II_list[0] * edge_num + (1 + L_list[1] + II_list[1]) * node_num
    print("Estimated exec time of agg_feat_in alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # update/apply phase for the aggrated results (phase 1)
    fifo_rst_agg_p1 = simpy.Store(env, capacity=16)
    L_list = [2, 13, 2]
    II_list = [1, 1, 1]
    env.process(update_agg(env, fifo_ft_h_agg, fifo_nod_src_list[3], fifo_rst_agg_p1, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (4 + L_list[0] + II_list[0] + L_list[1] + (64-1) * II_list[1] + L_list[2] + (16 - 1) * II_list[2]) * node_num
    print("Estimated exec time of update_agg alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))


    # update/apply phase for the aggrated results (phase 2)
    fifo_rst_agg = simpy.Store(env, capacity=10)
    L_list = [2, 9, 2]
    II_list = [1, 1, 1]
    env.process(update_agg_sum(env, fifo_rst_agg_p1, fifo_nod_src_list[4], fifo_rst_agg, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (4 + L_list[0] + (16 - 1) * II_list[0] + L_list[1] + (128-1)*II_list[1] + L_list[2] + II_list[2]) * node_num
    print("Estimated exec time of update_agg_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read the target feature from memory
    fifo_ft_in_tar = simpy.Store(env, capacity=10)
    L = 76
    II = 1
    L_mem = 64
    env.process(read_feat_in_tar(env, fifo_ft_in_tar, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_mem + L + II) * node_num
    print("Estimated exec time of read_feat_in_tar alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # update/apply phase for the target feature (phase 1)
    fifo_rst_tar_p1 = simpy.Store(env, capacity=16)
    L_list = [2, 13, 2]
    II_list = [1, 1, 1]
    env.process(update_tar(env, fifo_ft_in_tar, fifo_rst_tar_p1, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + II_list[0] + L_list[1] + (64-1) * II_list[1] + L_list[2] + (16 -1) * II_list[2]) * node_num
    print("Estimated exec time of update_tar alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # update/apply phase for the target feature (phase 2)
    fifo_rst_tar = simpy.Store(env, capacity=10)
    L_list = [2, 9, 2]
    II_list = [1, 1, 1]
    env.process(update_tar_sum(env, fifo_rst_tar_p1, fifo_rst_tar, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + (16 - 1) * II_list[0] + L_list[1] + (128-1)*II_list[1] + L_list[2] + II_list[2]) * node_num
    print("Estimated exec time of update_tar_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # concat the update/apply results of aggregated features and the target feature
    fifo_rst_cat = simpy.Store(env, capacity=10)
    L = 10
    II = 1
    env.process(concat_update_rst(env, fifo_rst_agg, fifo_rst_tar, fifo_nod_src_list[5], fifo_rst_cat, nidx_begin, nidx_end, L, II))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3 + L + II) * node_num
    print("Estimated exec time of concat_update_rst alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # normalization for each vertex (cal base + norm div)
    fifo_rst_cat2 = simpy.Store(env, capacity=10)
    fifo_rst_norm_reduce_sum = simpy.Store(env, capacity=16)
    L = 113
    II = 1
    env.process(norm_reduce_sum(env, fifo_rst_cat, fifo_rst_cat2, fifo_rst_norm_reduce_sum, nidx_begin, nidx_end, L, II))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + II * (node_num - 1)
    print("Estimated exec time of norm_reduce_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # accumulate norm reduced sum results of all FIFO segments
    fifo_rst_norm_acc = simpy.Store(env, capacity=10)
    L = 6
    II = 3
    env.process(norm_acc_rst(env, fifo_rst_norm_reduce_sum, fifo_rst_norm_acc, nidx_begin, nidx_end, L, II))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L + II + 10) * node_num
    print("Estimated exec time of norm_acc_rst alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # divide concat results by norm accumulate results
    fifo_rst_norm = simpy.Store(env, capacity=10)
    L = 18
    II = 1
    env.process(norm_div(env, fifo_rst_cat2, fifo_rst_norm_acc, fifo_rst_norm, nidx_begin, nidx_end, L, II))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + II * (node_num - 1)
    print("Estimated exec time of norm_div alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # write norm results to memory
    L_list = [2, 3]
    II_list = [1, 1]
    L_mem = 64
    env.process(write_norm_rst_mem(env, fifo_rst_norm, nidx_begin, nidx_end, L_list, II_list, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + II_list[0] + L_mem + L_list[1] + II_list[1]) * node_num
    print("Estimated exec time of write_norm_rst_mem alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))


    # run the simulation
    env.run()

    # measure the actural running time (seconds)
    end = timeit.default_timer()
    print("Actual running time: %s Seconds"%(end-start))

    # measure the simulation resultes (cycles)
    print("Simulation ends at %d" % env.now)
    print("Calculated simulation time (s): ", (env.now * 1000 / freq * 1e-9))


