#ifdef COOT_BUILD_WINDOWS

// Forcefully undefine some WinAPI defines that clash with symbols
// named the same way as these WinAPI macros.
//
// When including this file to work around namespace pollution caused
// by inclusion of Windows.h somewhere in the tree of files you include,
// mind that there are multiple problematic situations:
//
// Say you have a function "CollidingName()" and that Windows.h defines "CollidingName" as a macro
// that expands to "CollidingNameA".
// I.   If Windows.h gets included _after_ the declaration of "CollidingName()",
//      compilation will fail "No function CollidingNameA()" on the line where you
//      call "CollidingName()". The macro overwrote "CollidingName" to "CollidingNameA"
//      and there is no function "CollidingNameA()".
//
// II.  If Windows.h is included _before_ the declaration of "CollidingName()",
//      all occurences of "CollidingName" will be overwritten to "CollidingNameA"
//      and the code will compile correctly.
//      Mind that while the code will compile, you may run into linking issues
//      if the "CollidingName()" function comes from a library or a different translation unit.
//
// III. A rather tricky scenario:
//      1) Windows.h is included.
//      2) "CollidingName()" is declared.
//      3) This inhibitory header file gets included.
//      4) "CollidingName()" is caller or defined.
//      In this case, compilation will fail with "No function CollidingName()". In this case
//      the macro overwrote the declaration of "CollidingName()" to "CollidingNameA()" but we then
//      undefined the macro. Had we not done that, we would get the case II. behavior. If we
//      called "CollidingNameA()" instead, the code would compile, leaving everyone rather baffled.
//
// - If "CollidingName()" has split declaration and definition, the number of renaming issues
//   grows even larger.
// - These renaming issues apply to all of your code, class names, class members, templates, ..
//
// ***  Long story short - the order of includes, defines and undefs matters! ***

#undef GetAtomName
#undef IGNORE
#undef near
#undef far

// Tell others that we have inhibited parts of WinAPI
#ifndef COOT_WINAPI_INHIBITED
#define COOT_WINAPI_INHIBITED
#endif // COOT_WINAPI_INHIBITED

#endif // COOT_BUILD_WINDOWS
