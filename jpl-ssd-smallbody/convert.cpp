//g++ -shared -Llua convert.cpp -o convert.so
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

extern "C" {
int convert_int(lua_State* L) { return convert_T<int>(L); }
int convert_long(lua_State* L) { return convert_T<long>(L); }
int convert_float(lua_State* L) { return convert_T<float>(L); }
int convert_double(lua_State* L) { return convert_T<double>(L); }
}
