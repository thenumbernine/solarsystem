/*
compile with: g++ -shared -Llua convert.cpp -o convert.so
library for converting strings of raw binary data to numbers
*/
#include <lua.hpp>

template<typename T>
int lua_pushtype(lua_State* L, T value);

template<>
int lua_pushtype<int>(lua_State* L, int value) {
	lua_pushinteger(L, value);
	return 1;
}

template<>
int lua_pushtype<long>(lua_State* L, long value) {
	lua_pushinteger(L, value);
	return 1;
}

template<>
int lua_pushtype<float>(lua_State* L, float value) {
	lua_pushnumber(L, value);
	return 1;
}

template<>
int lua_pushtype<double>(lua_State* L, double value) {
	lua_pushnumber(L, value);
	return 1;
}

template<typename T>
int convert_T(lua_State* L) {
	const char* str = luaL_checkstring(L, 1);
	return lua_pushtype<T>(L, *(T*)str);
}

#ifdef REG
#error here
#endif
#define REG(x)	{#x, convert_T<x>}
static luaL_Reg funcs[] = {
	REG(int),
	REG(long),
	REG(float),
	REG(double),
};
#undef REG

extern "C" {
int luaopen_convert(lua_State* L) {
	luaL_newlib(L, funcs);
	return 1;
}
}
