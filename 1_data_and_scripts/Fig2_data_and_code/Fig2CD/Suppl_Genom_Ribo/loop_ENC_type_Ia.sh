#!/bin/bash
        for i in $( ls ); do
            echo item: $i
            chips -seqall $i -outfile ../1_ENC_Type_Ia/$i.chips
        done


