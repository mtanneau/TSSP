const TSSPDIR = "data/tssp"
const LOGDIR = "exp/log"

function print_job(fcore, ftime, fsto, frmp, flog)
    println("julia --project --sysimage=JuliaTSSP.so exp/colgen.jl $TSSPDIR/$fcore $TSSPDIR/$ftime $TSSPDIR/$fsto $frmp > $LOGDIR/$flog.log 2>&1")
    return nothing
end

# assets
print_job("assets.cor", "assets.tim", "assets.sto.large", "assets", "assets")

# 4node & 4node-base
for S in  2 .^ [10, 11, 12, 13, 14, 15]
    print_job("4node.cor", "4node.tim", "4node.sto.$S", "4node", "4node_$S")
    print_job("4node.cor.base", "4node.tim", "4node.sto.$S", "4node-base", "4node-base_$S")
end

# env & env-diss
for S in  [1200, 1875, 3780, 5292, 8232, 32928]
    print_job("env.cor", "env.tim", "env.sto.$S", "env", "env_$S")
    print_job("env.cor.diss", "env.tim", "env.sto.$S", "env-diss", "env-diss_$S")
end

# phone-large
print_job("phone.cor", "phone.tim", "phone.sto", "phone", "phone")

# stormG2
print_job("stormG2.cor", "stormG2.tim", "stormG2_1000.sto", "stormG2", "stormG2_1000")