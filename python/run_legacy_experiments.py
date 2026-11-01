"""Optional unnormalized integer baseline; not the intended bounded MAC."""
import random, statistics
from gamma_exact_model import encode_int, multiply_ints, dot_product_ints, span, nnz


def pct(vals,p):
    vals=sorted(vals)
    return vals[int(round((len(vals)-1)*p))]


def main(seed=0):
    rng=random.Random(seed)
    spans=[]; nz=[]
    for n in range(-128,128):
        d=encode_int(n); spans.append(span(d)); nz.append(nnz(d))
    print('UNNORMALIZED BASELINE: INT8 representation')
    print(' median span:', statistics.median(spans))
    print(' max span:', max(spans))
    print(' median nnz:', statistics.median(nz))

    psp=[]; tops=[]; rounds=[]; dprods=[]
    for _ in range(2000):
        a=rng.randint(-128,127); b=rng.randint(-128,127)
        result,st,raw,da,db=multiply_ints(a,b)
        psp.append(span(result)); tops.append(st.topplings); rounds.append(st.rounds); dprods.append(nnz(da)*nnz(db))
    print('\nINT8 products')
    print(' median span:', statistics.median(psp), 'p95:', pct(psp,.95))
    print(' median digit-products:', statistics.median(dprods))
    print(' median topplings:', statistics.median(tops), 'p95:', pct(tops,.95))
    print(' median rounds:', statistics.median(rounds), 'p95:', pct(rounds,.95))

    print('\nUNNORMALIZED BASELINE: 256-term dot normalization sweep (paired inputs)')
    workloads = [
        ([rng.randint(-128,127) for _ in range(256)],
         [rng.randint(-128,127) for _ in range(256)])
        for _ in range(20)
    ]
    for interval in [1,2,4,8,16,32,64,128,256]:
        tpm=[]; rpm=[]; peak=[]
        for xs, ws in workloads:
            result,st=dot_product_ints(xs,ws,normalize_every=interval)
            tpm.append(st.topplings_per_mac); rpm.append(st.rounds_per_mac); peak.append(st.peak_raw_coeff)
        print(f'every {interval:3d}: top/MAC={statistics.median(tpm):7.3f} rounds/MAC={statistics.median(rpm):7.4f} peak={statistics.median(peak):7.1f}')

if __name__=='__main__': main()
