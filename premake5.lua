require("Premake/Common")

workspace("MTL")
	common:setConfigsAndPlatforms()

	common:addCoreDefines()

	cppdialect("C++20")
	rtti("Off")
	exceptionhandling("Off")
	flags("MultiProcessorCompile")

	startproject("Tests")

	project("MTL")
		location("MTL/")
		warnings("Extra")

		kind("StaticLib")
		common:outDirs(true)

		common:addPCH("%{prj.location}/Src/PCH.cpp", "%{prj.location}/Src/PCH.h")

		includedirs({
			"%{prj.location}/Inc/",
			"%{prj.location}/Src/"
		})

		files({
			"%{prj.location}/Inc/**",
			"%{prj.location}/Src/**"
		})
		removefiles({ "*.DS_Store" })

		common:addActions()

	project("Tests")
		location("Tests/")
		warnings("Extra")

		common:outDirs()
		common:debugDir()

		filter("configurations:Debug")
			kind("ConsoleApp")

		filter("configurations:not Debug")
			kind("WindowedApp")

		filter({})

		includedirs({ "%{prj.location}/Src/" })

		files({ "%{prj.location}/Src/**" })
		removefiles({ "*.DS_Store" })

		links({ "MTL" })
		sysincludedirs({ "%{wks.location}/MTL/Inc/" })

		common:addActions()