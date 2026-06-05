# WHAMPT ICRH 仿真

## 安装

使用 Git 拉取仓库：

```cmd
git clone https://github.com/merplateau/VASIMR.git
```

## 推荐使用方式

### 准备配置

准备一份 `.nml` 配置文件，并将其所在文件夹路径填入 `config.json` 中。

```json
"nmlDir": "D:\\Data\\nml\\"
```

将 `compile.bat` 中的 Visual Studio 安装路径和 Intel OneAPI 环境变量启动脚本 `setvars.bat` 路径设置为实际路径：

```bat
set "VS2022INSTALLDIR=D:\Program Files\Microsoft Visual Studio\2022\Community"

call "D:\Program Files (x86)\Intel\oneAPI\setvars.bat"
```

### 确定版本号

将 `.nml` 文件命名为：`1.2.3@4.nml` 的形式，其中 @ 前的是代码版本号，@ 后的是该版本下的算例编号。

### 编译与运行

在项目根目录下，打开 `cmd`，运行：

```cmd
compile 1.2.3
```

编译中间产物会生成在 `\tmp\` 目录下，可执行文件会生成在 `\bin\` 目录下。

运行 case4：

```cmd
run 1.2.3@4
```

注意：使用 `run.bat` 快捷运行时只需要 `.nml` 文件民，会自动使用其名字中的代码版本号对应的可执行文件，请确保该版本已经被编译。

## 说明

本仓库可以手动或者使用 IDE 编译运行，但未经测试。

默认编译器为 ifx，可使用 ifort。