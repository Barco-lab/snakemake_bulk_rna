# bigwig：全部 reads 一个；可选按转录本正/负链拆分


rule bigwig:
    input:
        bam=final_bam,
    output:
        "results/07_bigwig/{sample}.bw",
    params:
        args=config["bigwig"]["params"],
    threads: threads_for("bigwig")
    log:
        "<logs>/bigwig/{sample}.log",
    benchmark:
        "<benchmarks>/bigwig/{sample}.tsv"
    conda:
        ENV_DEEPTOOLS
    shell:
        "bamCoverage -b {input.bam} -o {output} -p {threads} {params.args} > {log} 2>&1"


rule bigwig_stranded:
    """
    fwd/rev 指转录本所在链。deeptools 的 --filterRNAstrand 以 dUTP（reverse）文库为前提，
    forward 文库需要对调；unstranded 样本无法区分链，fwd/rev 都输出全部 reads。
    """
    input:
        bam=final_bam,
        s="results/04_count/strand/{sample}.s",
    output:
        fwd="results/07_bigwig/{sample}.fwd.bw",
        rev="results/07_bigwig/{sample}.rev.bw",
    params:
        args=config["bigwig"]["params"],
    threads: threads_for("bigwig")
    log:
        "<logs>/bigwig/{sample}.stranded.log",
    conda:
        ENV_DEEPTOOLS
    shell:
        """
        case $(cat {input.s}) in
            2) f="--filterRNAstrand forward"; r="--filterRNAstrand reverse" ;;
            1) f="--filterRNAstrand reverse"; r="--filterRNAstrand forward" ;;
            *) f=""; r=""; echo "WARNING: unstranded 文库，fwd/rev 均为全部 reads" > {log} ;;
        esac
        bamCoverage -b {input.bam} -o {output.fwd} -p {threads} {params.args} $f >> {log} 2>&1
        bamCoverage -b {input.bam} -o {output.rev} -p {threads} {params.args} $r >> {log} 2>&1
        """
