# Makefile for Fortran project

# コンパイラの指定
FC = gfortran

# コンパイルフラグ
FFLAGS = -Jbuild -Ibuild

# ソースファイル
SRCS = para.f90 calcu.f90 record.f90 main.f90

# オブジェクトファイル
OBJS = $(patsubst %.f90,build/%.o,$(SRCS))

# 実行ファイル名
EXEC = test

# デフォルトターゲット
all: build $(EXEC)

# 実行ファイルの生成ルール
$(EXEC): $(OBJS)
	$(FC) -o $@ $(OBJS)

# 各 .o ファイルの生成ルール
build/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# build ディレクトリが存在しない場合に作成するルール
build:
	mkdir -p build

# クリーンアップルール
clean:
	rm -f build/*.o build/*.mod $(EXEC)

# 実行ルール
run: $(EXEC)
	./$(EXEC)

.PHONY: all clean run build
