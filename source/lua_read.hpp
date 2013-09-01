#ifndef LUA_READ_H

#define LUA_READ_H

#include <stdio.h>
#include <string.h>

#include <iostream>
#include <map>
#include <lua.hpp>
//#include "/global/home/users/kigumen/modules/lua/5.2.1/include/lua.hpp"
//#include "/global/home/users/kigumen/modules/lua/5.1.5/include/lua.hpp"
//extern "C" { 
      ////#include <lua.h> 
      ////#include <lauxlib.h> 
      ////#include <lualib.h> 
      //#include "/global/home/users/kigumen/modules/lua/5.2.1/include/lua.h"
      ////#include "/global/home/users/kigumen/modules/lua/5.2.1/include/luaxlib.h"
      //#include "/global/home/users/kigumen/modules/lua/5.2.1/include/lualib.h"
//} 


using namespace std;


class Lua_read
{
private:
  lua_State *ls;
public:
  map<string,double> read_config(vector<string> vv) {
    map<string,double> ret;
    for(auto i:vv) {
      lua_getglobal(ls,i.c_str());
      ret.insert(make_pair(i,lua_tonumber(ls,-1)));
    }
    //for(auto i:vv) ret.insert(make_pair(i,lua_tonumber(i.c_str())));
    return ret;
  }
  Lua_read(const char *lua_config_filename) {
    ls = lua_open();
    luaL_openlibs(ls);
    if (luaL_loadfile(ls, lua_config_filename) || lua_pcall(ls, 0,0,0)) {
      cout << "error: " << lua_tostring(ls,-1) << "\n";
    }
    lua_pushnil(ls);
  }

  //template <class T>
  //T get(const char *var_name) {
    //char temp[64];
    //memset(temp, 0, sizeof(temp));
    //int i=0;
    //int j=0;
    //int n=0;
    //while (var_name[i] != '\0') {
      //char c = var_name[i];
      //if (c == '.') {
        //if (n == 0)
          //lua_getglobal(ls, temp);
        //else
          //lua_getfield(ls, -1, temp);
        //++n;
        //memset(temp, 0, sizeof(temp));
        //j = 0;

        //if (lua_isnil(ls, -1))
          //return 0;
      //}
      //else {
        //temp[j] = c;
        //++j;
      //}
      //++i;
    //}
    //if (n == 0)
      //lua_getglobal(ls, temp);
    //else
      //lua_getfield(ls, -1, temp);
    ////lua_getglobal(ls,var_name);
    //return lua_get<T>();
    //lua_pop(ls,1);
  //}

  //// Generic get
  //template<typename T>
  //T lua_get() {
    //return 0;
  //}
  // Specializations
};

  //template <> float         Lua_read::lua_get<float>() { return lua_tonumber(ls, -1); }
  //template <> double        Lua_read::lua_get<double>() { return lua_tonumber(ls, -1); }
  //template <> bool          Lua_read::lua_get<bool>() { return lua_toboolean(ls, -1); }
  //template <> int           Lua_read::lua_get<int>() { return lua_tointeger(ls, -1); }
  //template <> unsigned int  Lua_read::lua_get<unsigned int>() { return (unsigned int)lua_tonumber(ls, -1); }
  //template <> const char *  Lua_read::lua_get<const char *>() { return lua_tostring(ls, -1); }
  //template <> string Lua_read::lua_get<string>() { return string(lua_tostring(ls, -1)); }

#endif /* end of include guard: LUA_READ_H */

///// improve with boost sting http://www.boost.org/doc/libs/1_49_0/doc/html/string_algo.html

//////   Usage:
//
//#include <iostream>
//#include <lua.hpp>

//#include "lua_read.h"

//using namespace std;

//int main(int argc, char** argv)
//{
  //Lua_read lr("config.lua");

  //bool use_blur = lr.get<bool>("config.render.effects.blur");

  //cout << "vol is " << lr.get<double>("vol") << "\n";
  //cout << "rho is " << lr.get<int>("rho") << "\n";
  //cout << "method is " << lr.get<string>("method") << "\n";
  //cout << "method is " << lr.get<const char*>("method") << "\n";

  //return 0;
//}

/////>g++ main.C -I/usr/include/lua5.1/ -llua5.1 ; ./a.out                  
////http://stackoverflow.com/questions/7210154/unable-to-find-lua-headers-with-find-package-in-cmake
