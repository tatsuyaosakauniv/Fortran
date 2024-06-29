# Makefile for Fortran project

# コンパイラの指定
FC = gfortran

# コンパイルフラグ
FFLAGS = -Jbuild -Ibuild -Wall -Wextra -MMD -cpp

# ソースファイル
SRCS = para.f90 calcu.f90 record.f90 main.f90

# モジュールファイルの依存関係
MODS = $(patsubst %.f90,build/%.mod,$(SRCS))

# オブジェクトファイル
OBJS = $(patsubst %.f90,build/%.o,$(SRCS))

# 依存関係ファイル
DEPS = $(OBJS:.o=.d)

# 実行ファイル名
EXEC = test

# デフォルトターゲット
all: build $(EXEC)

# 実行ファイルの生成ルール
$(EXEC): $(OBJS)
	$(FC) -o $@ $(OBJS)

# 各 .o ファイルの生成ルール（依存関係ファイルのインクルード）
build/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
	@$(FC) $(FFLAGS) -MM -MT $@ $< > build/$*.d

# モジュールファイルの生成ルール
build/%.mod: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# build ディレクトリが存在しない場合に作成するルール
build:
	mkdir -p build

# クリーンアップルール
clean:
	rm -f build/*.o build/*.mod $(EXEC) $(DEPS)

# 実行ルール
run: $(EXEC)
	./$(EXEC)

# 依存関係ファイルのインクルード
-include $(DEPS)

.PHONY: all clean run build
