#pragma once
#include <map>
#include <string>
#include <chrono>
#include <iostream>
#include <stack>

class MyTimer
{
public:
	typedef std::chrono::duration<double> t_type;
	MyTimer(const std::string &name) :name(name), isRunning(true){
		start();
	};
//	MyTimer(const char *name) :name(name), isRunning(true){
//		start();
//	};
	~MyTimer() {
		end();
	};
	void start() {
		startTime = std::chrono::system_clock::now();
		isRunning = true;
	}
	void end() {
		if(isRunning){
			endTime = std::chrono::system_clock::now();
			t_type t = endTime - startTime;
			// printf(" %s: %f\n", name.c_str(), t.count());
			if(name != ""){
				_tmMap()[name] += t;
				_tmMapCnt()[name]++;
			}
			isRunning = false;
		}
	}
	double getTime(){
		return (endTime-startTime).count();
	} 

	template<typename F, typename... Args>
	static t_type funcTime(F&& func, Args&&... args){
		auto t0 = std::chrono::system_clock::now();
		func(std::forward<Args>(args)...);
		auto t1 = std::chrono::system_clock::now();
		return t_type(t1-t0);
	}

	static void pusht() {
		auto t = std::chrono::system_clock::now();
		return get_tstack().push(t);
	}
	static double popt() {
		auto t = get_tstack().top();
		get_tstack().pop();
        auto duration = std::chrono::duration_cast<t_type>(std::chrono::system_clock::now() - t);
		return duration.count();
	}

	static void clear() {
		_tmMap().clear();
		_tmMapCnt().clear();
	}
	static void clear(const std::string &name) {
		_tmMap()[name] = t_type::zero();
	}
	static double get(const std::string &name) {
		return _tmMap()[name].count();
	}

	static std::stack< std::chrono::time_point<std::chrono::system_clock> >& get_tstack()
	{
		static std::stack< std::chrono::time_point<std::chrono::system_clock> > _tstack;
		return _tstack;
	}
	

//	static double get(const char *name) {
//		//std::string n(name);
//		return g_tmMap[std::string(name)].count();
//	}
//	static int cnt(const char *name) {
//		return g_tmMapCnt[std::string(name)];
//	}
	static int cnt(const std::string &name) {
		return _tmMapCnt()[std::string(name)];
	}
//	static void print(const char *name, std::ostream &outp = std::cout) {
//		outp << "Time for " << name << ": " << get(name) << "  avg:" << get(name) / cnt(name) << std::endl;
//	}
	static void print(const std::string &name, std::ostream &outp = std::cout) {
		outp << "Time for " << name << ": " << get(name) << "  avg:" << get(name) / cnt(name) << "  cnt:" << cnt(name) << std::endl;
	}
	static void printAll(std::ostream &outp = std::cout) {
		for (auto it = _tmMap().begin(); it != _tmMap().end(); it++) {
			print(it->first, outp);
		}
		clear();
	}
	static void printCurTime(std::ostream &outp = std::cout){
		auto time_point = std::chrono::system_clock::now();
	    std::time_t ttp = std::chrono::system_clock::to_time_t(time_point);
		outp << ttp << std::endl;
	}

    template<typename F, typename ...Args>
    static double measure(F func, Args&&... args) {
        auto start = std::chrono::system_clock::now();

        func(std::forward<Args>(args)...);

        auto duration = std::chrono::duration_cast<t_type>(std::chrono::system_clock::now() - start);

        return duration.count();
    }

private:
	static std::map<std::string, t_type> &_tmMap(){
		static std::map<std::string, t_type> g_tmMap;
		return g_tmMap;
	}
	static std::map<std::string, int> &_tmMapCnt(){
		static std::map<std::string, int> g_tmMapCnt;
		return g_tmMapCnt;
	}
//	static std::map<std::string, t_type> g_tmMap;
//	static std::map<std::string, int> g_tmMapCnt;
	std::chrono::time_point<std::chrono::system_clock> startTime;
	std::chrono::time_point<std::chrono::system_clock> endTime;
	std::string name;
	bool isRunning;

};
