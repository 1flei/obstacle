#pragma once
#include <set>
#include <unordered_map>
#include "subtidx.h"

//priority queue used for subtidx

class SubtPriQueue
{
    using KeyType = std::pair<double, SubtIdx>;
    std::set<KeyType> key_set;
    std::unordered_map<int, std::set<KeyType>::iterator > iter_map;

public:
    void emplace(double dist, SubtIdx st) {
        int tid = st.tid;
        auto it = iter_map.find(tid);
        if(it != iter_map.end()){
            //tid is already in queue
            if(dist < it->second->first){
                key_set.erase(it->second);
                iter_map[tid] = key_set.emplace(dist, st).first;
            }
        } else{
            iter_map[tid] = key_set.emplace(dist, st).first;
        }
    }

    void emplace_k(double dist, SubtIdx st, int k) {
        if(size() < k) {
            emplace(dist, st);
        } else{
            auto& maxk = maximum();
            if(dist < maxk.first){ 
                //push, will either update the current key or push a new key
                emplace(dist, st);
                if(size() > k){
                    //if pushed a new key
                    pop_maximum();
                }
            }
        }
    }

    size_t size(){
        return key_set.size();
    }

    const KeyType& top() {
        return *(key_set.begin());
    }

    const KeyType& minimum() {
        return *(key_set.begin());
    }

    const KeyType& maximum() {
        return *(key_set.rbegin());
    }

    void pop() {
        auto kt = top();
        iter_map.erase(kt.second.tid);
        key_set.erase(key_set.begin());
    }

    void pop_minimum() {
        pop();
    }

    void pop_maximum() {
        auto kt = maximum();
        iter_map.erase(kt.second.tid);
        key_set.erase(std::prev(key_set.end()));
    }

    std::set<KeyType>::iterator begin() {
        return key_set.begin();
    }
    std::set<KeyType>::iterator end() {
        return key_set.end();
    }
    const std::set<KeyType>::iterator begin() const {
        return key_set.begin();
    }
    const std::set<KeyType>::iterator end() const {
        return key_set.end();
    }

    bool empty() const {
        return key_set.empty();
    }

    void erase(const KeyType& kt){
        key_set.erase(kt);
        iter_map.erase(kt.second.tid);
    }

    std::set<KeyType>::iterator erase(std::set<KeyType>::iterator it){
        auto ret = key_set.erase(it);
        iter_map.erase(it->second.tid);
        return ret;
    }
};