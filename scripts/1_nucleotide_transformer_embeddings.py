from transformers import AutoTokenizer, AutoModelForMaskedLM
import torch
from Bio import SeqIO
import time
import numpy as np

#### INPUTS ####
enhancers_fasta = ''
output = 'data/NT_embeddings.csv'
#####

device = torch.device("cuda")
cpu = torch.device('cpu')

tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True)
model = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True)
model = model.to(device)

fasta_sequences = SeqIO.parse(open(enhancers_fasta),'fasta')

all_embeddings = []

t = time.time()
for i,x in enumerate(fasta_sequences):
    seq = str(x.seq)
    if i % 1000 == 0:
        print(f'Time elapsed for 1000 sequences: {time.time() - t:.1f} seconds.')
        t = time.time()
    tokens_ids = tokenizer.batch_encode_plus([seq], return_tensors="pt", padding="max_length", max_length = 351).input_ids
    tokens_ids = tokens_ids.to(device)
    attention_mask = tokens_ids != tokenizer.pad_token_id
    torch_outs = model(
            tokens_ids,
            attention_mask=attention_mask,
            encoder_attention_mask=attention_mask,
            output_hidden_states=True,
    )
    embeddings = torch_outs['hidden_states'][-1].detach()
    attention_mask = torch.unsqueeze(attention_mask, dim=-1)
    mean_sequence_embeddings = (torch.sum(attention_mask*embeddings, axis=-2)/torch.sum(attention_mask, axis=1)).flatten().to(cpu)
    all_embeddings.append(mean_sequence_embeddings)

np.savetxt(output, torch.stack(all_embeddings).numpy(), fmt='%.5f', delimiter=',')
