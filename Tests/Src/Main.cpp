#include <MTL/Utils/Core.h>
#include <MTL/ComplexNumber.h>

int main(int argc, char** argv)
{
	constexpr MTL::ComplexNumber complex = MTL::SqrtC(-4.0f);

}

#if BUILD_IS_SYSTEM_WINDOWS

#include <Windows.h>

int WinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ LPSTR lpszCmdLine, _In_ int nCmdShow)
{
	return main(__argc, __argv);
}

#endif