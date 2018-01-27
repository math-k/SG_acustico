#!/bin/bash

#Variavel que guarda inicio do tempo de execucao
start_time=$(date +%s)
clear

#Espacamento entre tiros, e receptores
spacing="$(sed '12q;d' input.txt)"
spacing="$(cut -d' ' -f1 <<<"$spacing")"
#Posicao original da fonte lida de input.txt
source_position="$(sed '4q;d' input.txt)"
source_position="$(cut -d' ' -f1 <<<"$source_position")"
og_sp=$source_position
#Posicao original do primeiro receptor do streamer lida de input.txt
streamer_tail="$(sed '10q;d' input.txt)"
streamer_tail="$(cut -d' ' -f1 <<<"$streamer_tail")"
og_st=$streamer_tail
#Posicao original do ultimo receptor do streamer lida de input.txt
streamer_head="$(sed '11q;d' input.txt)"
streamer_head="$(cut -d' ' -f1 <<<"$streamer_head")"
og_sh=$streamer_head
#Numero de tiros lido de input.txt
n_shots="$(sed '9q;d' input.txt)"
n_shots="$(cut -d' ' -f1 <<<"$n_shots")"
echo "Numero de tiros=$n_shots"
#Comprimento do streamer
streamer_length="$(sed '13q;d' input.txt)"
streamer_length="$(cut -d' ' -f1 <<<"$streamer_length")"
echo "Comprimento do streamer=$streamer_length"

#Variavel que guarda arquivo de input do script stack.f90
INPUT_STACK="/home/mathk/Desktop/git/SG_acustico/Stack/input_stack.txt"

#Grava o numero de tiros(numero de imagens) no arquivo de input do script de stack
echo $n_shots >> $INPUT_STACK

#variaveis auxiliares que guardam a posicao anterior da fonte, do primeiro e do ultimo receptor do streamer
aux=0
aux2=0
aux3=0
#Loop que varia a posicao da fonte e do streamer
for ((i=1; i<=n_shots; i++)); do
    echo "Tiro: $i" 
    ./onda_2d    
    aux=$source_position
    aux2=$streamer_tail
    aux3=$streamer_head
    source_position=$((source_position + spacing))
    streamer_tail=$((streamer_tail + spacing))
    streamer_head=$((streamer_head + spacing))
    sed -i "4s/$aux/$source_position/" input.txt
    sed -i "10s/$aux2/$streamer_tail/" input.txt
    sed -i "11s/$aux3/$streamer_head/" input.txt
    mv seismogram.bin seismogram$i.bin
    mv image.bin image$i.bin
done
#Copia as imagens geradas para o diretorio do script de stack
cp image*.bin Stack
#Remove as imagens do diretorio atual
rm image*.bin 
#Muda para o diretorio do script de stack
cd Stack
#Empilhamento das imagens
./stack


#Retorna os arquivos de input para o estado inicial e retorna para o diretorio inicial
sed -i '2d' $INPUT_STACK
cd ..
sed -i "4s/$source_position/$og_sp/" input.txt
sed -i "10s/$streamer_tail/$og_st/" input.txt
sed -i "11s/$streamer_head/$og_sh/" input.txt

#Calculo do tempo de execucao
end_time=$(date +%s)
run_time=$(((end_time-start_time)/60))
echo "Runtime: $run_time minutes."
