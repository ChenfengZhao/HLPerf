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

def read_pseudo_in_mem(env, fifo_in_nod_src, fifo_out_pseudo_in, nidx_begin, nidx_end, L, II, L_mem):
    """read pseudo_in_mat from memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_pseudo_in : simpy.store
        Output FIFO of pseudo of each edge
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
            
            yield fifo_out_pseudo_in.put(1)

    print("read_pseudo_in_mem ends at %d" % env.now)

def update_pseudo_pp(env, fifo_in_pseudo_in, fifo_in_nod_src, fifo_out_pseudo_pp, nidx_begin, nidx_end, L, II):
    """compute pseudo_pp by pseudo_proj (matrix multiplication)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_pseudo_in : simpy.store
        Input FIFO of pseudo of each edge
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_pseudo_pp : simpy.store
        Output FIFO of generated pseudo of each edge
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
            _ = yield fifo_in_pseudo_in.get()
            
            if e == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)

            yield fifo_out_pseudo_pp.put(1)
    
    print("update_pseudo_pp ends at %d" % env.now)

def comp_weight_gaussian(env, fifo_in_pseudo_pp, fifo_in_nod_src, fifo_out_gaussian, nidx_begin, nidx_end, L, II):
    """compute gaussian weight

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_pseudo_pp : simpy.store
        Input FIFO of generated pseudo of each edge
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_gaussian : simpy.store
        Output FIFO of gaussian results of each edge
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
            _ = yield fifo_in_pseudo_pp.get()
            
            if e == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            yield fifo_out_gaussian.put(1)

    print("comp_weight_gaussian ends at %d" % env.now)

def read_feat_in_agg(env, fifo_in_nod_src, fifo_in_tmp_src, fifo_out_ft_in_agg, nidx_begin, nidx_end, L, II, L_mem=64):
    """read input features to be aggregated from memory

    Instead of putting acutal input feature of nbrs, we put feke values

    Parameters  
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_in_tmp_src :simpy.store
        Input FIFO of neighbor indices of each node
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
    L_mem : int
        the latency of memory read access
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        # yield env.timeout(L_mem) # memory latency
        # yield env.timeout(15) # memory latency

        for e in range(deg):
            _ = yield fifo_in_tmp_src.get()

            if e == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)

            yield fifo_out_ft_in_agg.put(1)
                
    
    print("read_feat_in_agg ends at %d" % env.now)

def agg_feat_in(env, fifo_in_gaussian, fifo_in_nod_src, fifo_in_ft_in_agg, fifo_out_rst_agg, nidx_begin, nidx_end, L, II):
    """aggregate features while multiply with e_gaussian 
    ft_in_1 * e_gaussian_1 + ft_in_2 * e_gaussian_2 + ....

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_gaussian : simpy.store
        Input FIFO of gaussian results of each edge
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_in_ft_in_agg : simpy.store
        Input FIFO of input feature to be aggretated
    fifo_out_rst_agg : simpy.store
        Output FIFO of aggretaed fature
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
        yield env.timeout(4)

        # if deg != 0:
        for e in range(deg):
            _ = yield fifo_in_ft_in_agg.get()
            _ = yield fifo_in_gaussian.get()

            if e == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            if e == (deg - 1):
                yield fifo_out_rst_agg.put(1)
                
        
            # yield env.timeout(40) # estimated latency since the loops are unrolled
            # yield fifo_out_rst_agg.put(1)

    print("agg_feat_in ends at %d" % env.now)


def update_agg_fc(env, fifo_in_rst_agg, fifo_in_nod_src, fifo_out_rst_update_p1, nidx_begin, nidx_end, L_list, II_list):
    """update aggregated results using fc weights (phase 1)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_agg : simpy.store
        Input FIFO of aggretaed fature
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_update_p1 : simpy.store
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
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(4)

        if deg != 0:
            _ = yield fifo_in_rst_agg.get()
            yield env.timeout(1)

            for i in range(2):
                if i == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
            
            yield env.timeout(L_list[1] + (64-1) * II_list[1])

            for kd in range(8):
                if kd == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
                
                yield fifo_out_rst_update_p1.put(1)

def update_agg_fc_sum(env, fifo_in_rst_update_p1, fifo_in_nod_src, fifo_out_rst_update, nidx_begin, nidx_end, L_list, II_list):
    """update aggregated results using fc weights (phase 1)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_update_p1 : simpy.store
        Input FIFO of result feature (phase 1) of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_update : simpy.store
        Output FIFO of result feature (phase 2) of each node
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
            for kd in range(8):
                _ = yield fifo_in_rst_update_p1.get()
                if kd == 0:
                    yield env.timeout(L_list[0])
                else:
                    yield env.timeout(II_list[0])
                
            yield env.timeout(L_list[1] + (64-1)*II_list[1])

            for i in range(2):
                if i == 0:
                    yield env.timeout(L_list[2])
                else:
                    yield env.timeout(II_list[2])
                
                yield fifo_out_rst_update.put(1)
        
    print("update_agg_fc_sum ends at %d" % env.now)

def ksum_rst_update(env, fifo_in_rst_update, fifo_in_nod_src, fifo_out_rst_ksum, nidx_begin, nidx_end, L, II):
    """sum rst_update in the order of kernel

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_update : simpy.store
        Input FIFO of result feature (phase 2) of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_ksum : simpy.store
        Output FIFO of result ksum feature of each node
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
        yield env.timeout(4)

        if deg != 0:
            for i in range(2):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)
                
                _ = yield fifo_in_rst_update.get()
            
            yield env.timeout(5)

            yield fifo_out_rst_ksum.put(1)


def write_ksum_rst_mem(env, fifo_in_rst_ksum, fifo_in_nod_src, nidx_begin, nidx_end, L, II, L_mem=64):
    """write ksum results to memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_ksum : simpy.store
        Input FIFO of result ksum feature of each node
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
            _ = yield fifo_in_rst_ksum.get()
            yield env.timeout(L_mem) # memory latency
            if n == nidx_begin:
                yield env.timeout(L)
            else:
                yield env.timeout(II)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    add_option(parser)
    args = parser.parse_args()

    # read input graph
    format = "gnnhls"

    # input_fn = "./data/OGBG-MOLTOX21/csr_indptr_trans.txt"
    # input_fn = "./data/OGBG-MOLHIV/csr_indptr_trans.txt"
    # input_fn = "./data/OGBN-ARXIV/csr_indptr_trans.txt"
    # input_fn = "./data/OGBN-PROTEINS/csr_indptr_trans.txt"

    # input_fn = input_fn_dict[args.ds_fn]
    input_fn = gen_input_fn(args.ds_fn)
    print("Input dataset:", args.ds_fn, ":", input_fn)

    # hardware parameters (MHz)
    freq = 300
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
    fifo_nod_src_list = [simpy.Store(env, capacity=10), simpy.Store(env, capacity=11), simpy.Store(env, capacity=12), simpy.Store(env, capacity=13), simpy.Store(env, capacity=14), simpy.Store(env, capacity=15), simpy.Store(env, capacity=16), simpy.Store(env, capacity=17), simpy.Store(env, capacity=18), simpy.Store(env, capacity=19)]
    L = 4
    II = 1
    L_mem = 64
    env.process(read_nod_src(env, fifo_nod_src_list, nidx_begin, nidx_end, L, II, L_mem, deg_list))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L + L_mem + II * (node_num - 1)
    print("Estimated exec time of read_nod_src alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read edge_src from memory
    fifo_tmp_src = simpy.Store(env, capacity=14)
    L = 75
    II = 1
    L_mem = 64
    env.process(read_edge_src(env, fifo_nod_src_list[0], fifo_tmp_src, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + L_mem * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of read_edge_src alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read pseudo_in_mat from memory
    fifo_pseudo_in = simpy.Store(env, capacity=10)
    L = 75
    II = 1
    L_mem = 64
    env.process(read_pseudo_in_mem(env, fifo_nod_src_list[1], fifo_pseudo_in, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + L_mem * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of read_pseudo_in_mem alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))


    # compute pseudo_pp by pseudo_proj (matrix multiplication)
    fifo_pseudo_pp = simpy.Store(env, capacity=10)
    L = 103
    II = 1
    env.process(update_pseudo_pp(env, fifo_pseudo_in, fifo_nod_src_list[2], fifo_pseudo_pp, nidx_begin, nidx_end, L, II))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of update_pseudo_pp alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # compute gaussian weight
    fifo_gaussian = simpy.Store(env, capacity=11)
    L = 69
    II = 1
    env.process(comp_weight_gaussian(env, fifo_pseudo_pp, fifo_nod_src_list[3], fifo_gaussian, nidx_begin, nidx_end, L, II))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of comp_weight_gaussian alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read input features to be aggregated from memory
    fifo_ft_in_agg = simpy.Store(env, capacity=10)
    L = 77
    II = 2
    L_mem = 64
    env.process(read_feat_in_agg(env, fifo_nod_src_list[4], fifo_tmp_src, fifo_ft_in_agg, nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3 + L_mem) * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of read_feat_in_agg alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # aggregate features while multiply with e_gaussian (ft_in_1 * e_gaussian_1 + ft_in_2 * e_gaussian_2 + ....)
    fifo_rst_agg = simpy.Store(env, capacity=10)
    L = 16
    II = 4
    env.process(agg_feat_in(env, fifo_gaussian, fifo_nod_src_list[5], fifo_ft_in_agg, fifo_rst_agg, nidx_begin, nidx_end, L, II))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of agg_feat_in alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # update aggregated results
    fifo_rst_update_p1 = simpy.Store(env, capacity=32)
    L_list = [1, 12, 2]
    II_list = [1, 1, 1]

    env.process(update_agg_fc(env,fifo_rst_agg, fifo_nod_src_list[6], fifo_rst_update_p1, nidx_begin, nidx_end, L_list, II_list))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + 1 * node_num + (L_list[0] + II_list[0] + L_list[1] + (64-1) * II_list[1] + L_list[2] + (8-1) * II_list[2]) * node_num
    print("Estimated exec time of update_agg_fc alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # update aggregated results using fc weights (phase 1)
    fifo_rst_update = simpy.Store(env, capacity=10)
    L_list = [2, 9, 2]
    II_list = [1, 1, 1]

    env.process(update_agg_fc_sum(env, fifo_rst_update_p1, fifo_nod_src_list[7], fifo_rst_update, nidx_begin, nidx_end, L_list, II_list))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + (L_list[0] + (8-1) * II_list[0] + L_list[1] + (64-1)*II_list[1] + L_list[2] + II_list[2]) * node_num
    print("Estimated exec time of update_agg_fc_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # sum rst_update in the order of kernel
    fifo_rst_ksum = simpy.Store(env, capacity=10)
    L = 10
    II = 4
    env.process(ksum_rst_update(env, fifo_rst_update, fifo_nod_src_list[8], fifo_rst_ksum, nidx_begin, nidx_end, L, II))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + (L + II + 5) * node_num
    print("Estimated exec time of ksum_rst_update alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # write ksum results to memory
    L = 74
    II = 2
    L_mem = 64
    env.process(write_ksum_rst_mem(env, fifo_rst_ksum, fifo_nod_src_list[9], nidx_begin, nidx_end, L, II, L_mem))

    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 4 * node_num + L_mem * node_num + L + II * (node_num - 1)
    print("Estimated exec time of write_ksum_rst_mem alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # run the simulation
    env.run()

    # measure the actural running time (seconds)
    end = timeit.default_timer()
    print("Actual running time: %s Seconds"%(end-start))

    # measure the simulation resultes (cycles)
    print("Simulation ends at %d" % env.now)
    print("Calculated simulation time (s): ", (env.now * 1000 / freq * 1e-9))
    






