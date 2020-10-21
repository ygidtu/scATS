#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
后修改的Bed类，满足自己的需求.

2017.3.6 进一步进行修改和完善，补充了正负链的判断过程
2017.4.20 完成了对bed类的修改
共计两个Bed类
共计两个Bed类

Bed，赋予Bed文件的路径，及自动读取，并且可以进行相应的计算
（1）包含的内置参数
        __read_bed: 自动调用，读取Bed文件，强制sort
        __save_bed：根据指定的类型，输出bed文件，
                    output: 指定输出文件，如果不指定，默认输出到stdout
                    contain: 指定输出某个beds，有[contain，hits，not_hits]
        __save_both: 根据指定的类型，同时输出两个对应的bed
                     output: 同上
                     contain: 只有[contain, hits]
        __save_closest: 保存closest的结果
                        output：同上
        __get_index: 读取完bed文件后，按照染色体生成字典，bed6是sort过的, {chr: [Bed6]}
        __reverse_dict: 针对获得的hits和contain来将key， value反过来的，用来赋值给第二个Bed

        参数说明：
            contain: 包含关系的数据，contain, c
            hits: 具有重合位点的数据， hits, hit, h
            not_hits: 不具有任何重合位点的数据，not_hits， not_hit， not, n
（2）重载的神奇方法
        __add__: 将两个Bed接起来，并且sort
        __sub__: 取第一个Bed独有的Bed
        __and__：找Bed6之间具有重合位点的，生成词典，key是该Bed的，value是第二个Bed的
        __or__：找Bed6之间完全不具有重合位点的，list格式，Bed6
        __len__: 返回Bed的长度，即list的长度
        __contains__: 判断Bed6是否包含在该Bed中
        __reversed__: 将Bed中bed的顺序反过来
        __hash__：利用copy将bed复制了一下，防止顺序被打乱，sort之后，提取所有的bed，计算hash
（3）添加的自定义方法
        cover: Bed1中完全覆盖Bed2的，返回字典
        closest: Bed1和Bed2中位点的距离，center（默认False），返回字典{Bed1: [[Bed2, dist]]}
        save: 保存文件用的
                output: 输出文件，没有，就重新输出到stdout
                save: 要保存的数据类型 bed（单个bed），both（两个bed同时），closest（距离信息）
                contain: 同上
        distance: 计算两个bed之间距离，
                  center（默认False），用法同上
                  dist，只保留该距离内的数据，默认为4000
（4）内置数据
        beds: 读取玩bed之后，由Bed6组成的list
        hits: 具有重合位点的Bed的信息，{Bed1: [Bed2]}
        contain: Bed1完全覆盖Bed2的信息， {Bed1: [Bed2]}
        not_hits: 没有任何重合位点的部分bed， list
        dist：计算过后的两个Bed之间的距离 {Bed1: [[Bed2, dist]]}

"""
import sys
import os
from copy import deepcopy
from utils.bed6 import Bed6


class Bed(object):
    u"""Bed类，真正供外界使用的类."""

    def __init__(self, bed=None):
        u"""初始化."""
        if bed is None:
            self.beds = []
        elif isinstance(bed, list):
            if all([isinstance(x, Bed6) for x in bed]):
                self.beds = bed
            else:
                self.beds = []
                for bed_ in bed:
                    if os.path.exists(bed_):
                        self.beds += self.__read_bed(bed)
        else:
            self.beds = self.__read_bed(bed)
        self.beds.sort()

        self.__hits = self.beds
        self.__contain = None
        self.__not_hits = None
        self.__dist = None

    def add(self, bed: Bed6):
        self.beds.append(bed)

    def __read_bed(self, bed):
        u"""读取bed文件的."""

        header = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'append']
        beds = []
        with open(bed) as reader:
            for line in reader:
                if line.startswith('chr'):
                    # 正常打开文件需要进行一个判断
                    # 根据每行内容的多少，生成相应的bed
                    lines = line.split()
                    # 读取的6列基本构成一个完整的字典文件
                    tem = {key: value for key,
                           value in zip(header[:6], lines[:6])}
                    # 添加完整的附录信息
                    tem.update({header[-1]: lines[6:]})

                    beds.append(self.dict_to_bed(tem))

        # beds.sort()
        return beds

    def dict_to_bed(self, indict):
        u"""将字典转化成Bed6."""
        return Bed6(chrom=indict.get('chrom'), start=indict.get('start'),
                    end=indict.get('end'), name=indict.get('name'),
                    score=indict.get('score'), strand=indict.get('strand'),
                    append=indict.get('append'))

    '''
    保存文件或者输出系列
    '''

    def sort(self):
        u"""手动进行sort()吧，加减之后无法sort."""
        self.beds = sorted(self.beds)
        return self

    def self_merge(self, strandness=True):
        u"""
        合并自身重复区域
        :param strandness: 是否严格根据不同链进行修正
        """
        if len(self.beds) < 2:
            return self

        data = {}
        for i in self.beds:
            key = i.chrom if strandness else "%s#%s" % (i.chrom, i.strand)
            temp = data.get(key, [])
            temp.append(i)
            data[key] = temp

        res = []
        for beds in data.values():
            last = beds[0]

            for i in range(1, len(beds)):
                if last > i:
                    res.append(last)
                    last = i
                else:
                    last += i
            
            res.append(last)
        
        res.sort()
        self.beds = res

    def reduce(self):
        u"""减少除去重复项."""
        self.beds = list(set(self.beds))
        self.__hits = deepcopy(self.beds)
        return self

    '''
    保存系列function
    '''

    def __save_bed(self, hits, output=None):
        u"""将bed文件输出的."""
        if output is None:
            for i in hits:
                print(i.get_bed())
        else:
            with open(output, 'w+') as writer:
                for i in hits:
                    writer.write(i.get_bed() + '\n')

    def __save_both(self, hits, output=None):
        u"""将bed文件输出的，同时输出两个文件的."""
        if output is None:
            for key, values in hits.items():
                for value in values:
                    sys.stdout.write(format('%s|%s\n' % (key.get_bed(),
                                                         value.get_bed())))
        else:
            with open(output, 'w+') as w:
                for key, values in hits.items():
                    for value in values:
                        w.write(format('%s|%s\n' % (key.get_bed(),
                                                    value.get_bed())))

    def save_cover(self, output=None, both=False):
        u"""用来保存cover结果的."""
        if self.__contain:
            if both:
                self.__save_both(self.__contain, output=output)
            else:
                self.__save_bed(self.__contain, output=output)

    def save_closest(self, output=None):
        u"""用来保存closest结果的."""

        def construct_(values):
            u"""构成一行结果，按照|分割."""
            return '|'.join(['{0}|{1}'.format(value.get_bed(),
                                              dist) for value, dist in values])

        u"""
        2017.8.18
        Python缝缝补补又三年的作用终于发挥出来了
        downstream_target由于需要特殊的判断正负链，
        因此此处，特地拆分成正负链两个最近，以满足需求

        但是，没有分两个链的部分，会分别在两个链的结果中出现，
        这样子就无法使用原本的字典来储存结果，因此改用list

        然而，由于使用了同一套输出系统，而这里为了适应list的变化，
        特地添加了一个判断过程，分别针对这两个做处理
        """
        if output is not None and self.__dist:
            if not isinstance(self.__dist, list):
                dist = [self.__dist]
            else:
                dist = self.__dist
            with open(output, 'w+') as w:
                for dist_ in dist:
                    for this, value in dist_.items():
                        w.write(format('%s|%s\n' %
                                       (this.get_bed(), construct_(value))))
        elif self.__dist:
            if not isinstance(self.__dist, list):
                dist = [self.__dist]
            else:
                dist = self.__dist

            for dist_ in dist:
                for this, value in dist_.items():
                    sys.stdout.write(
                        format(
                            '%s|%s\n' %
                            (
                                this.get_bed(), construct_(value)
                            )
                        )
                    )

    def save(self, output=None, both=False, hits=True):
        u"""用来保存文件的函数，可以选择保存bed或者both."""
        if hits and self.__hits:
            if both:
                self.__save_both(self.__hits, output=output)
            else:
                self.__save_bed(self.__hits, output=output)
        else:
            hits = self.__not_hits
            self.__save_bed(hits, output=output)

    '''
    写几个getters，方便取值
    '''

    def get_hits(self):
        u"""返回hits的结果."""
        return self.__hits

    def get_not_hits(self):
        u"""返回没有hits的结果."""
        return self.__not_hits

    def get_cover(self):
        u"""返回完全覆盖的位点的信息."""
        return self.__contain

    def get_closest(self):
        u"""返回距离最近的所有数据."""
        return self.__dist

    '''
    2017.7.3

    返回所有的bed的某一样属性，返回一个list列表
    '''

    def get_beds(self):
        u"""返回所有的bed数据."""
        return [x.get_bed() for x in self.beds]

    def get_names(self):
        u"""返回所有bed的名字，方便用来进行对比确定都有啥."""
        return [x.get_name() for x in self.beds]

    def get_chroms(self):
        u"""返回所有的bed的染色体."""
        return [x.get_chrom() for x in self.beds]

    def get_starts(self):
        u"""返回所有的bed的起始位点."""
        return [x.get_start() for x in self.beds]

    def get_ends(self):
        u"""返回所有bed的终止位点."""
        return [x.get_end() for x in self.beds]

    def get_strands(self):
        u"""返回所有的bed的链信息."""
        return [x.get_strand() for x in self.beds]

    '''
    我的算法不在基于index来进行分类，因此直接将index删除即可
    '''

    def __reverse_dict(self, dict_):
        u"""将dict的key和value反过来，value默认是list."""
        hits = {}
        for key, values in dict_.items():
            tem = set([key])
            for value in values:
                if value in hits.keys():
                    tem = tem | hits[value]
                hits.update({value: tem})
        # 将原本的set转换成list
        for key, value in hits.items():
            value = list(value)
            hits.update({key: value})
        return hits

    def __add__(self, other):
        u"""将两个bed合并到一起的."""
        self.beds.extend(other.beds)
        self.beds.sort()
        return self.beds

    def __sub__(self, other):
        u"""去除第二个中含有的bed，只保留该bed中独有的."""
        self.beds = list(set(self.beds) - set(other.beds))
        self.beds.sort()
        return self.beds

    def __and__(self, other):
        u"""重载&, 用来找有一个重合位点的."""
        self.__hits = {}

        i, j = 0, 0

        while i < len(self.beds) and j < len(other.beds):
            # 算法重写，如果该文件的位点在另一个文件的上游，那么就将该文件的指针往后挪一个
            if self.beds[i] < other.beds[j]:
                i += 1
            # 匹配中了，就保留结果，然后该文件的指针后移一位
            elif self.beds[i] & other.beds[j]:
                tem = [other.beds[j]]
                if self.beds[i] in self.__hits.keys():
                    tem.extend(self.__hits[self.beds[i]])

                self.__hits.update({self.beds[i]: tem})
                i += 1
            # 其他的所有异常情况，就全部是另一个问价的指针往后挪一位
            else:
                j += 1

        '''
        反向生成other的hits，但是这个结果并不准确，无法包含所有值
        没有什么好的解决方案，重新跑一边的话，何苦呢.

        保留的唯一意义可能就是用于找完全没有任何结合位点的
        '''
        other.__hits = self.__reverse_dict(self.__hits)
        return self.__hits

    def __or__(self, other):
        u"""用来获取，完全没有任何重合的位点的."""
        '''
        调用内置函数，然后先找到有重合位点的，随后将其除掉，保留剩余的
        调用方法与普通的方法一致
        '''
        self.__and__(other)
        self.__not_hits = list(set(self.beds) - set(self.__hits.keys()))
        other.__not_hits = list(set(other.beds) - set(other.__hits.keys()))
        return self.__not_hits

    def __len__(self):
        u"""返回自定义类的长度."""
        return len(self.beds)

    def __contains__(self, item):
        u"""in，判断某个Bed是否是in."""
        if item in self.beds:
            return True
        else:
            return False

    def __reversed__(self):
        u"""将beds的排列顺序取反."""
        self.beds.reverse()

    def __iter__(self):
        u"""生成器，便于在列表外部调用所有的数据来处理文件."""
        for bed in self.beds:
            yield bed

    def __hash__(self):
        u"""算一个hash值."""
        # 复制一下，省的将原本的顺序打乱掉
        from copy import copy
        tem = copy(self.beds)
        tem.sort()
        msg = ''
        for i in tem:
            msg += i.get_bed()
        return hash(msg)

    def cover(self, other):
        u"""乘法，用来找出在范围内的，完全包含关系的，本bed的是否完全覆盖另外的."""
        '''
        强制调用一次&，省的导入了不同的bed进行比对，然后出问题.
        这个也没啥好改进得了，反正原本的基础上，也是这么个算法
        '''
        self.__and__(other)
        self.__contain = {}
        # 能够进行覆盖度查询的，怎么也得是至少有一个重合位点的才行吧？
        for key, values in self.__hits.items():
            for value in values:
                if key * value:
                    tem = [value]
                    if key in self.__contain.keys():
                        tem.extend(self.__contain[key])
                    self.__contain.update({key: tem})
        # 对于other则是被包含的内容
        other.contain = self.__reverse_dict(self.__contain)
        return self.__contain

    def closest(self, other, closest_=True, center_=False):
        u"""
        类似bedops的closest-feature，
        只能输出start site(或center)之间距离最近的那一组
        （只输出没有任何重合位点的）.
        """
        up = self.__closest_up(other, center_)
        down, and_ = self.__closest_down(other, center_)

        def merge_two_dicts(fir, sec):
            u"""两个字典整合成一个."""
            common = sorted(list(set(fir.keys()) | set(sec.keys())))

            merged = {}
            for key in common:
                tem = []
                if key in fir.keys():
                    tem += fir[key]
                if key in sec.keys():
                    tem += sec[key]
                tem = sorted(tem, key=lambda x: x[1])
                merged.update({key: tem})
            return merged

        self.__dist = merge_two_dicts(up, merge_two_dicts(down, and_))

        if closest_:
            self.select_closest()
        return self.__dist

    def __closest_down(self, other, center_=False):
        u"""找到self的上游离这个脚本最近的那个位点."""
        down = {}
        and_ = {}       # 其中有重合位点的，距离为0

        i, j = 0, 0
        while i < len(self.beds) and j < len(other.beds):
            statu = self.beds[i].is_same_chrom(other.beds[j])

            if statu is True:
                # 计算默认的距离
                distance = other.beds[j].get_end() - self.beds[i].get_start()
                # 如果是center，就计算区域中心的距离
                if center_:
                    distance = other.beds[j] / self.beds[i]
                # 如果该文件的bed处于另外一个bed的下游，那么该文件的指针往后挪一位
                if self.beds[i] > other.beds[j]:
                    down.update(
                        {self.beds[i]: [(other.beds[j], distance - 1)]})
                    j += 1
                # 如果两个文件具有重合位点了，就将另外一个文件往后推一位
                elif self.beds[i] & other.beds[j]:
                    and_.update({self.beds[i]: [(other.beds[j], 0)]})
                    i += 1
                # 除此之外的情况，就另一个文件往后推一位，这个另外情况，只有该文件在另一个文件的下游这种情况了吧
                else:
                    i += 1
            elif statu == 1:
                j += 1
            elif statu == -1:
                i += 1

        return down, and_

    def __closest_up(self, other, center_=False):
        u"""
        找到self的下游最近的那个点.
        算法与找上游的一样，区别就在于反向便利
        """
        dist = {}
        i, j = len(self.beds) - 1, len(other.beds) - 1
        while i >= 0 and j >= 0:

            statu = self.beds[i].is_same_chrom(other.beds[j])

            if statu is True:
                distance = other.beds[j].get_start() - self.beds[i].get_end()
                if center_:
                    distance = other.beds[j] / self.beds[i]
                if self.beds[i] < (other.beds[j]):
                    dist.update(
                        {self.beds[i]: [(other.beds[j], distance + 1)]})
                    j -= 1
                else:
                    i -= 1
            elif statu == 1:
                i -= 1
            elif statu == -1:
                j -= 1
        return dist

    def select_closest(self):
        u"""选择最近的结果."""
        for key, values in self.__dist.items():
            sort = sorted(values, key=lambda x: abs(x[1]))
            self.__dist.update({key: sort[:1]})

    def _split_pos_neg_(self):
        u"""将正负链分开."""
        pos = []
        neg = []

        for x in self.beds:
            if x.get_strand() == '-':
                neg.append(x)
            elif x.get_strand() == '+':
                pos.append(x)
            else:
                neg.append(x)
                pos.append(x)

        return sorted(pos), sorted(neg)

    def _downstream_real_(self, this, that, pos=True):
        u"""
        downstream_target正链部分的处理
        """
        dist = {}
        i, j = 0, 0
        while i < len(this) and j < len(that):
            # 判断两个区域的上下游关系
            statu = this[i].isdownstream_same(that[j])

            if statu is True:
                # 计算距离
                distance = int(that[j].get_score()) - this[i].get_center()

                # 如果是负链，负链越大的越上游，减出来的距离应当是正数
                if pos is False and distance >= 0:
                    dist.update({this[i]: [that[j], distance]})
                # 正链的话，数字越小的越上游，减出来的应当为负数
                elif pos is True and distance <= 0:
                    dist.update({this[i]: [that[j], distance]})

                i += 1
            elif statu == 2:
                j += 1
            elif statu == 3:
                i += 1
        return dist

    def downstream_target(self, other):
        u"""
        找下游最近的靶基因和他们的距离.
        要求：
        self的中心区域要在other的TSS位点的上游
        """
        # 提取两者共有的所有的bed，并且排序
        this_pos, this_neg = self._split_pos_neg_()
        that_pos, that_neg = other._split_pos_neg_()

        # 合并算法的插入算法，确实速度刚刚的
        results = [self._downstream_real_(this_pos, that_pos, pos=True)]

        results.append(self._downstream_real_(this_neg, that_neg, pos=False))

        return self.__dist


if __name__ == '__main__':
    pass
