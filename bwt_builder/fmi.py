from collections import defaultdict

from .bwt import build_bwt
from .bwt import build_suffix_array

class FMIndex:
    def __init__(self, text, sentinel="$", checkpoint_step=16):
        if sentinel not in text:
            text += sentinel
        self.sentinel = sentinel
        self.text = text
        self.sa = build_suffix_array(text)
        self.bwt = build_bwt(text, sentinel=sentinel)
        self.checkpoint_step = checkpoint_step
        self.first_occurrence = self.build_first_occurrence()
        self.checkpoints = self.build_checkpoints()

    def build_checkpoints(self):
        counts = defaultdict(list)
        total_counts = defaultdict(int)
        chars = set(self.bwt)
        for c in chars:
            counts[c] = []
        for i, char in enumerate(self.bwt):
            if i % self.checkpoint_step == 0:
                for c in chars:
                    counts[c].append((i, total_counts[c]))
            total_counts[char] += 1
        for c in chars:
            counts[c].append((len(self.bwt), total_counts[c]))
        return counts

    def build_first_occurrence(self):
        first_col = sorted(self.bwt)
        first_occurrence = {}
        for i, c in enumerate(first_col):
            if c not in first_occurrence:
                first_occurrence[c] = i
        return first_occurrence
    
    def count_symbol(self, symbol, pos):
        if symbol not in self.checkpoints:
            return 0
        if pos <= 0:
            return 0
        pos = min(pos, len(self.bwt))
        
        checkpoints = self.checkpoints[symbol]
        step = self.checkpoint_step
        checkpoint_idx = pos // step
        
        if checkpoint_idx < len(checkpoints):
            checkpoint_pos, checkpoint_count = checkpoints[checkpoint_idx]
        else:
            checkpoint_pos, checkpoint_count = checkpoints[-1]
        count = checkpoint_count
        
        for i in range(checkpoint_pos, pos):
            if self.bwt[i] == symbol:
                count += 1
        
        return count
    
    def search_exact(self, pattern):
        if not pattern:
            return list(range(len(self.text) - 1))
        
        top = 0
        bottom = len(self.bwt) - 1

        for symbol in reversed(pattern):
            if symbol not in self.first_occurrence:
                return []
            
            top = self.first_occurrence[symbol] + self.count_symbol(symbol, top)
            bottom = self.first_occurrence[symbol] + self.count_symbol(symbol, bottom + 1) - 1
            
            if top > bottom:
                return []
        
        result = []
        for i in range(top, bottom + 1):
            sa_val = self.sa[i]
            if sa_val < len(self.text) - 1:
                result.append(sa_val)
        
        return result